#include <iostream>
#include <vector>
#include <string>
#include <thread>
#include <chrono>
#include <cstring>
#include <png.h>

void grayscaleSequential(png_bytep* image, int width, int height, int channels) {
    for (int y = 0; y < height; ++y) {
        png_bytep row = image[y];
        for (int x = 0; x < width; ++x) {
            png_bytep px = &(row[x * channels]);
            uint8_t gray = static_cast<uint8_t>(0.3 * px[0] + 0.59 * px[1] + 0.11 * px[2]);
            px[0] = px[1] = px[2] = gray;
        }
    }
}

//final exam - implement this function
extern void grayscaleThreaded(png_bytep* image, int width, int height, int channels, int numThreads) ;

png_bytep* deepCopyImage(png_bytep* image, int height, int rowbytes) {
    png_bytep* copy = new png_bytep[height];
    for (int y = 0; y < height; ++y) {
        copy[y] = new png_byte[rowbytes];
        std::memcpy(copy[y], image[y], rowbytes);
    }
    return copy;
}

void freeImage(png_bytep* image, int height) {
    for (int y = 0; y < height; ++y)
        delete[] image[y];
    delete[] image;
}

void writeImage(const char* baseOutputFile, const char* suffix, png_bytep* image, int width, int height) {
    std::string outputPath(baseOutputFile);
    size_t dotPos = outputPath.rfind('.');
    std::string finalPath = (dotPos != std::string::npos) ?
        outputPath.substr(0, dotPos) + suffix + outputPath.substr(dotPos) :
        outputPath + suffix + ".png";

    FILE* output = fopen(finalPath.c_str(), "wb");
    if (!output) {
        std::cerr << "Error opening output file " << finalPath << "\n";
        return;
    }

    png_structp pngWrite = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    png_infop infoWrite = png_create_info_struct(pngWrite);
    if (!pngWrite || !infoWrite || setjmp(png_jmpbuf(pngWrite))) {
        std::cerr << "Error initializing PNG write for " << finalPath << "\n";
        if (pngWrite) png_destroy_write_struct(&pngWrite, &infoWrite);
        fclose(output);
        return;
    }

    png_init_io(pngWrite, output);
    png_set_IHDR(pngWrite, infoWrite, width, height, 8,
                 PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_write_info(pngWrite, infoWrite);
    png_write_image(pngWrite, image);
    png_write_end(pngWrite, nullptr);

    png_destroy_write_struct(&pngWrite, &infoWrite);
    fclose(output);
    std::cout << "Wrote output: " << finalPath << "\n";
}

int main(int argc, char* argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " input.png output.png [numThreads]\n";
        return 1;
    }

    const char* inputFile = argv[1];
    const char* outputFile = argv[2];
    int numThreads = (argc > 3) ? std::stoi(argv[3]) : std::thread::hardware_concurrency();

    FILE* fp = fopen(inputFile, "rb");
    if (!fp) {
        std::cerr << "Error opening input file\n";
        return 1;
    }

    png_structp png = png_create_read_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (!png) return 1;

    png_infop info = png_create_info_struct(png);
    if (!info) return 1;

    if (setjmp(png_jmpbuf(png))) return 1;

    png_init_io(png, fp);
    png_read_info(png, info);

    int width = png_get_image_width(png, info);
    int height = png_get_image_height(png, info);
    png_byte color_type = png_get_color_type(png, info);
    png_byte bit_depth = png_get_bit_depth(png, info);

    if (bit_depth == 16)
        png_set_strip_16(png);
    if (color_type == PNG_COLOR_TYPE_PALETTE)
        png_set_palette_to_rgb(png);
    if (color_type == PNG_COLOR_TYPE_GRAY && bit_depth < 8)
        png_set_expand_gray_1_2_4_to_8(png);

    if (png_get_valid(png, info, PNG_INFO_tRNS))
    png_set_tRNS_to_alpha(png);
    if (color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA) png_set_gray_to_rgb(png);
    if (color_type == PNG_COLOR_TYPE_RGB_ALPHA || color_type == PNG_COLOR_TYPE_GRAY_ALPHA || png_get_valid(png, info, PNG_INFO_tRNS)) 
    {
    	// Strip alpha channel if present
    	png_set_strip_alpha(png);
    }
    png_read_update_info(png, info);

    int channels = png_get_channels(png, info);
    if (channels != 3) {
    std::cerr << "Expected RGB image with 3 channels, but got " << channels << "\n";
    return 1;
}

    png_bytep* image = new png_bytep[height];
    for (int y = 0; y < height; ++y)
        image[y] = new png_byte[png_get_rowbytes(png, info)];

    png_read_image(png, image);
    fclose(fp);

    size_t rowbytes = png_get_rowbytes(png, info);
    png_bytep* imageSeq = deepCopyImage(image, height, rowbytes);
    png_bytep* imageThreaded = deepCopyImage(image, height, rowbytes);
    freeImage(image, height);

    // Time sequential grayscale
    auto t0 = std::chrono::high_resolution_clock::now();
    //_grayscaleSequential(imageSeq, width, height);
    grayscaleSequential(imageSeq, width, height, channels);
    auto t1 = std::chrono::high_resolution_clock::now();
    double seqTime = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "Sequential grayscale time: " << seqTime << " ms\n";

    // Time threaded grayscale
    auto t2 = std::chrono::high_resolution_clock::now();
    //_grayscaleThreaded(imageThreaded, width, height, numThreads);
    grayscaleThreaded(imageThreaded, width, height, channels, numThreads);
    auto t3 = std::chrono::high_resolution_clock::now();
    double threadTime = std::chrono::duration<double, std::milli>(t3 - t2).count();
    std::cout << "Threaded grayscale time (" << numThreads << " threads): " << threadTime << " ms\n";

    // Write output images
    writeImage(outputFile, "_seq", imageSeq, width, height);
    writeImage(outputFile, "_thread", imageThreaded, width, height);

    freeImage(imageSeq, height);
    freeImage(imageThreaded, height);
    return 0;
}
// kenroche@amd.com
// k8r@uw.edu
