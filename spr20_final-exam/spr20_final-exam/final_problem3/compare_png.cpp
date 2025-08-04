#include <iostream>
#include <vector>
#include <string>
#include <png.h>
#include <cstdlib>

bool loadPNG(const std::string& filename, png_bytep*& row_pointers, int& width, int& height) {
    FILE* fp = fopen(filename.c_str(), "rb");
    if (!fp) {
        std::cerr << "Failed to open file: " << filename << "\n";
        return false;
    }

    png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (!png) return false;

    png_infop info = png_create_info_struct(png);
    if (!info) return false;

    if (setjmp(png_jmpbuf(png))) return false;

    png_init_io(png, fp);
    png_read_info(png, info);

    width  = png_get_image_width(png, info);
    height = png_get_image_height(png, info);
    png_byte color_type = png_get_color_type(png, info);
    png_byte bit_depth  = png_get_bit_depth(png, info);

    if (bit_depth == 16) png_set_strip_16(png);
    if (color_type == PNG_COLOR_TYPE_PALETTE) png_set_palette_to_rgb(png);
    if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8) png_set_expand_gray_1_2_4_to_8(png);
    if (png_get_valid(png, info, PNG_INFO_tRNS)) png_set_tRNS_to_alpha(png);
    if (color_type == PNG_COLOR_TYPE_RGB || color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_PALETTE)
        png_set_filler(png, 0xFF, PNG_FILLER_AFTER);
    if (color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA)
        png_set_gray_to_rgb(png);

    png_read_update_info(png, info);

    row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * height);
    for (int y = 0; y < height; y++) {
        row_pointers[y] = (png_byte*)malloc(png_get_rowbytes(png, info));
    }

    png_read_image(png, row_pointers);
    fclose(fp);
    png_destroy_read_struct(&png, &info, nullptr);
    return true;
}

void freePNG(png_bytep* row_pointers, int height) {
    for (int y = 0; y < height; y++) free(row_pointers[y]);
    free(row_pointers);
}

bool comparePNGs(png_bytep* img1, png_bytep* img2, int width, int height) {
    bool identical = true;
    for (int y = 0; y < height; y++) {
        png_bytep row1 = img1[y];
        png_bytep row2 = img2[y];
        for (int x = 0; x < width * 4; x++) { // RGBA = 4 bytes/pixel
            if (row1[x] != row2[x]) {
                std::cout << "Difference at (y=" << y << ", x=" << x/4 << "), channel=" << (x % 4)
                          << ": " << (int)row1[x] << " vs " << (int)row2[x] << "\n";
                identical = false;
            }
        }
    }
    return identical;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: ./compare_png file1.png file2.png\n";
        return 1;
    }

    png_bytep *img1 = nullptr, *img2 = nullptr;
    int w1, h1, w2, h2;

    if (!loadPNG(argv[1], img1, w1, h1) || !loadPNG(argv[2], img2, w2, h2)) {
        std::cerr << "Error loading images.\n";
        return 1;
    }

    if (w1 != w2 || h1 != h2) {
        std::cerr << "Image dimensions do not match.\n";
        freePNG(img1, h1);
        freePNG(img2, h2);
        return 1;
    }

    bool ok = comparePNGs(img1, img2, w1, h1);
    std::cout << (ok ? "✅ Images are identical.\n" : "❌ Images differ.\n");

    freePNG(img1, h1);
    freePNG(img2, h2);
    return 0;
}

