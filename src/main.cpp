#include <cmath>
#include <iostream>
#include <vector>
#include <algorithm>
#include "lodepng.h"

using namespace std;

typedef unsigned char channel; // using unsigned chars to store 8-bit image channels in line with the lodepng standard

class seam_carver {
    public:
        seam_carver (char* filepath) // decode the input image into a vector<channel> and set the values of width and height
        {
            lodepng::decode(img, width, height, filepath);
        }

        void carve (int seams) // carve out one seam at a time in a loop
        {
            for (int i {0}; i < seams; ++i, --width)
                img = get_carved_image();
        }

        void output_image (char* filepath)
        {
            if (lodepng::encode(filepath, img, width, height)) {
                cerr << "Error: output image could not be decoded.\n";
                exit(2);
            }
        }

    private:
        vector<channel> img;
        unsigned int width, height;

        vector<channel> get_greyscale() // converts a vector of RGBA channels to greyscale pixels of a single channel
        {
            size_t size {width * height};
            vector<channel> ret (size);

            for (size_t i {0}; i < size; ++i) {
                size_t times_4 {4 * i};
                ret[i] = (img[times_4] + img[times_4 + 1] + img[times_4 + 2]) / 3; // since the original RGBA image has 4 channels and the greyscale image will have just one channel, the index of the ith pixel in the original image and its channels is 4 times the greyscale image's pixel index
            }

            return ret;
        }

        vector<vector<int>> get_sobel() // takes a vector of greyscale pixels and returns a 2D vector of ints containing the Sobel operator results for each pixel
        {
            auto img_greyscale {get_greyscale()};
            vector<vector<int>> ret (height, vector<int> (width));

            for (size_t y {0}; y < height; ++y) {
                bool has_up   {y != 0},
                     has_down {y != height - 1};

                for (size_t x {0}; x < width; ++x) {
                    bool has_left  {x != 0},
                         has_right {x != width - 1};

                    auto current_pos {y * width + x};

                    // NOTE: it's more efficient to use img_greyscale[current_pos] in these ternary operators instead of initialising an lvalue because all the pixels inside of the image border - which typically account for most of the image - will not need to use this variable
                    auto top_left     {has_left  ? (has_up   ? img_greyscale[current_pos - width - 1] : img_greyscale[current_pos - 1]) : img_greyscale[current_pos]},
                         top_right    {has_right ? (has_up   ? img_greyscale[current_pos - width + 1] : img_greyscale[current_pos + 1]) : img_greyscale[current_pos]},
                         bottom_left  {has_left  ? (has_down ? img_greyscale[current_pos + width - 1] : img_greyscale[current_pos - 1]) : img_greyscale[current_pos]},
                         bottom_right {has_right ? (has_down ? img_greyscale[current_pos + width + 1] : img_greyscale[current_pos + 1]) : img_greyscale[current_pos]};

                    // calculate the x image of the Sobel operator
                    int sobel_x {top_left - top_right + // UB is 255 * 4 = 1020
                                2 * ((has_left ? img_greyscale[current_pos - 1] : img_greyscale[current_pos]) - (has_right ? img_greyscale[current_pos + 1] : img_greyscale[current_pos])) +
                                bottom_left - bottom_right};

                    // calculate the y image of the Sobel operator
                    int sobel_y {top_left - bottom_left + // UB is also 255 * 4 = 1020
                                2 * ((has_up ? img_greyscale[current_pos - width] : img_greyscale[current_pos]) - (has_down ? img_greyscale[current_pos + width] : img_greyscale[current_pos])) +
                                top_right - bottom_right};

                    ret[y][x] = sqrt (sobel_x * sobel_x + sobel_y * sobel_y); // NOTE: the upper bound of the Sobel operator on a 24-bit pixel is 1442 (which is approx. sqrt (2*1020*1020))
                }
            }

            return ret;
        }

        vector<size_t> get_seam() // takes a 2D vector of ints and finds the lowest energy path in it via dynamic programming, using the Sobel operator output vector to store its previous values
        {
            auto dp {get_sobel()};

            size_t y_max {height - 1}, x_max {width - 1};
            vector<vector<size_t>> predecessor (y_max, vector<size_t> (x_max + 1)); // stores the column position of the lowest-cost ongoing path which the current coordinate came from, e.g. predecessor[0][0] contains the column position that dp[1][0] came from

            size_t final_min; // stores the column of the pixel in the final row with the minimum path value

            {
                size_t y {1};

                auto min_predecessor = [&dp, &y](size_t a, size_t b) { // comparison function that returns the index with a lesser energy path in the previous row. used to calculate the energy of the current pixel but also to mark its predecessor
                    return dp[y-1][a] < dp[y-1][b];
                };

                for (; y < y_max; ++y) { // the first row of pixels shouldn't change
                    auto current_min {min (0, 1, min_predecessor)};
                    dp[y][0] += dp[y-1][current_min];
                    predecessor[y-1][0] = current_min;

                    for (size_t x {1}; x < x_max; ++x) {
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

                for (size_t x {1}; x < x_max; ++x) {
                    {
                        auto min_pred {min ({x-1, x, x+1}, min_predecessor)};
                        dp[y][x] += dp[y-1][min_pred], // all the neighbours below and directly to the right pixel aren't checked because their minimum paths haven't been calculated yet
                        predecessor[y - 1][x] =  min_pred;
                    }

                    if (dp[y][final_min] > dp[y][x])
                        final_min = x;
                }

                if (dp[y][x_max] + min (dp[y-1][x_max - 1], dp[y-1][x_max]) < dp[y][final_min])
                    final_min = x_max;
            }

            vector<size_t> ret (y_max + 1); // contains the pixels that led to the minimum path
            ret[y_max] = final_min;

            for (size_t y {y_max - 1};; --y) {
                ret[y] = predecessor[y][final_min],
                final_min = ret[y];

                if (y == 0) // prevent bit overflow by checking if y == 0 before decrementing y again
                    return ret;
            }
        }

        vector<channel> get_carved_image() // outputs a vector<channel> for lodepng to encode, containing a version of the image with one seam carved out
        {
            auto seam {get_seam()};
            size_t row_channels {4 * width};
            vector<channel> ret (height * (row_channels - 4));

            for (size_t y {0}, img_row {0}; y < height; ++y, img_row += row_channels) // img_row represents the current row in the original image
                for (size_t x {0}, move_forward {0}, current_seam {4 * seam[y]};; ++x) {
                    if (x == current_seam) // if the current x value is part of the seam we should remove...
                        move_forward = 4; // then change the move_forward counter to 4 so we can move 1 RGBA pixel ahead and skip the pixel in the seam

                    if (x + move_forward >= row_channels) // we should check if we need to break before updating the image but also after updating the move_forward variable. if we don't do this then a heap corruption will occur when the program attempts to remove a seam that ends at the last pixel of a row
                        break;

                    ret[img_row - 4 * y + x] = img[img_row + x + move_forward]; // we take away 4y from img_row because we have skipped 4 channels by removing 1 pixel from each of the y preceding rows
                }

            return ret;
        }
};

int main(int argc, char** argv)
{
    if (argc != 4) {
        cerr << "Error: bad number of inputs. Please only input an input .png filename, an output .png filename, and a number of seams to remove.\n";
        return 1;
    }

    seam_carver sc (argv[1]); // input filepath
    sc.carve(atoi(argv[3])); // number of seams to carve out
    sc.output_image(argv[2]); // output filepath
}
