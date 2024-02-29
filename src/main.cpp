#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include "lodepng.h"

using namespace std;

// using unsigned chars to store integers in the range 0-255 in line with the lodepng standard
typedef unsigned char channel;

// converts a vector of RGBA channels to greyscale pixels with single channels
vector<channel> to_greyscale (vector<channel>& img, unsigned& width, unsigned& height) {
    vector<channel> ret (height * width);

    for (size_t y = 0, row_channels = 4 * width; y < height; ++y) {
        size_t rows_passed = y * row_channels;
        for (size_t x = 0; x < width; ++x) {
            auto current_pixel = rows_passed + 4 * x;
            ret[current_pixel / 4] = (img[current_pixel] + img[current_pixel + 1] + img[current_pixel + 2]) / 3;
        }
    }

    return ret;
}

// takes a vector of greyscale pixels and returns a 2D vector of ints containing the Sobel operator results for each pixel
// NOTE: although the upper bound of the Sobel operator fits snugly in an uint16_t, a 32-bit int is used instead in this function. this is because the get_seam function will have to use a 32-bit integer 2D vector of the same dimensions anyway as part of its dynamic programming approach to finding the lowest energy path. thanks to copy elision we can prevent having to initialise two different 2D vectors and save up on memory in the process
vector<vector<int>> to_sobel (vector<channel> img, unsigned& width, unsigned& height) {
    vector<vector<int>> ret (height, vector<int> (width));

    for (size_t y = 0; y < height; ++y) {
        bool has_up   = y != 0,
             has_down = y != height - 1;

        for (size_t x = 0; x < width; ++x) {
            bool has_left  = x != 0,
                 has_right = x != width - 1;

            auto current_pos = y * width + x;

            // NOTE: it's actually more efficient to use img[current_pos] in these ternary operators instead of initialising an lvalue because all the pixels inside of the image border - which typically account for most of the image - will not need to use this variable
            auto top_left     = has_left ? (has_up ? img[current_pos - width - 1] : img[current_pos - 1]) : img[current_pos],
                 top_right    = has_right ? (has_up ? img[current_pos - width + 1] : img[current_pos + 1]) : img[current_pos],
                 bottom_left  = has_left ? (has_down ? img[current_pos + width - 1] : img[current_pos - 1]) : img[current_pos],
                 bottom_right = has_right ? (has_down ? img[current_pos + width + 1] : img[current_pos + 1]) : img[current_pos];

            // calculates the x image of the Sobel operator
            int sobel_x = top_left - top_right + // UB is 255 * 4 = 1020
                          2 * ((has_left ? img[current_pos - 1] : img[current_pos]) - (has_right ? img[current_pos + 1] : img[current_pos])) +
                          bottom_left - bottom_right;

            // calculates the y image of the Sobel operator
            int sobel_y = top_left - bottom_left + // UB is also 255 * 4 = 1020
                          2 * ((has_up ? img[current_pos - width] : img[current_pos]) - (has_down ? img[current_pos + width] : img[current_pos])) +
                          top_right - bottom_right;

            ret[y][x] = sqrt (sobel_x * sobel_x + sobel_y * sobel_y); // NOTE: the upper bound of the Sobel operator on a 24-bit pixel is 1442 (which is approx. sqrt (2*1020*1020))
        }
    }

    return ret;
}

// takes a 2D vector of ints and finds the lowest energy path from it via dynamic programming, using the input 2D vector to store previous values
vector<size_t> get_seam (vector<vector<int>> dp) {
    size_t y_max = dp.size() - 1, x_max = dp[0].size() - 1;
    vector<vector<size_t>> predecessor (y_max, vector<size_t> (x_max + 1)); // stores the column position of the lowest-cost ongoing path which the current coordinate came from, e.g. predecessor[0][0] contains the column position that dp[1][0] came from

    size_t final_min; // stores the column of the pixel in the final row with the minimum path value

    {
        size_t y = 1;

        auto min_predecessor = [&dp, &y](size_t a, size_t b) { // comparison function that returns the index with a lesser energy path in the previous row. used to calculate the energy of the current pixel but also to mark its predecessor
            return dp[y-1][a] < dp[y-1][b];
        };

        for (; y < y_max; ++y) { // the first row of pixels shouldn't change
            auto current_min = min (0, 1, min_predecessor);
            dp[y][0] += dp[y-1][current_min];
            predecessor[y-1][0] = current_min;

            for (size_t x = 1; x < x_max; ++x) {
                current_min = min ({x-1, x, x+1}, min_predecessor);
                dp[y][x] += dp[y-1][current_min], // all the neighbours below and directly to the right pixel aren't checked because their minimum paths haven't been calculated yet. the pixel immediately to the left isn't checked either because this would also reduce the height of the image
                predecessor[y - 1][x] = current_min;
            }

            current_min = min (x_max - 1, x_max, min_predecessor);
            dp[y][x_max] += dp[y-1][current_min];
            predecessor[y-1][x_max] = current_min;
        }

        final_min = min (0, 1, min_predecessor);
        dp[y][0] += dp[y-1][final_min];

        for (size_t x = 1; x < x_max; ++x) {
            auto min_pred = min ({x-1, x, x+1}, min_predecessor);
            dp[y][x] += dp[y-1][min_pred], // all the neighbours below and directly to the right pixel aren't checked because their minimum paths haven't been calculated yet
            predecessor[y - 1][x] =  min_pred;

            if (dp[y][final_min] > dp[y][x])
                final_min = x;
        }

        if (dp[y][x_max] + min (dp[y-1][x_max - 1], dp[y-1][x_max]) < dp[y][final_min])
            final_min = std::move(x_max);
    }

    vector<size_t> ret (y_max + 1); // contains the pixels that led to the minimum path
    ret[y_max] = final_min;

    for (size_t y = y_max - 1;; --y) {
        ret[y] = predecessor[y][final_min],
        final_min = ret[y];

        if (y == 0) // prevent bit overflow by checking if y == 0 before decrementing y again
            return ret;
    }
}

// takes an image, a vector of size_ts represents the pixel taken by the minimum energy path at each row, and the images' dimensions to output a final, seam carved image
vector<channel> carve_seam (vector<channel>& img, vector<size_t> seam, unsigned& width, unsigned& height) {
    size_t row_channels = 4 * width,
           new_row_channels = row_channels - 4;

    vector<channel> ret (height * new_row_channels);

    for (size_t y = 0, ret_row = 0, img_row = 0; y < height; ++y, ret_row += new_row_channels, img_row += row_channels) // ret_row represents the current row in the return image, while img_row represents the current row in the original image. they differ in size because the return image has one less pixel than the original image in each of its rows.
        for (size_t x = 0, add = 0, current_seam = 4 * seam[y];; ++x) {
            if (x == current_seam) // if the current x value is part of the seam we should remove...
                add = 4; // then change the add counter to 4 so we can move 1 RGBA pixel ahead and skip the pixel in the seam

            if (x + add >= row_channels) // we should check if we need to break before updating the image but also after updating the "add" variable. if we don't do this then a heap corruption will occur when the program attempts to remove a seam that ends at a pixel at the end of its row
                break;

            ret[ret_row + x] = std::move(img[img_row + x + add]);
        }

    return ret;
}

int main(int argc, char** argv) {
    if (argc != 4) {
        cerr << "Error: bad number of inputs. Please only input an input .png filename, an output .png filename, and a number of seams to remove.\n";
        return 1;
    }

    vector<channel> img;
    unsigned width, height;
    if (lodepng::decode(img, width, height, argv[1])) {
        cerr << "Error: image could not be decoded. Have you entered a valid .png file?\n";
        return 2;
    }

    for (int i = 0, times = atoi(argv[3]); i < times; ++i, --width)
        img = carve_seam(img, get_seam (to_sobel(to_greyscale(img, width, height), width, height)), width, height);

    if (lodepng::encode(argv[2], img, width, height)) {
        cerr << "Error: output image could not be decoded.\n";
        return 3;
    }
}
