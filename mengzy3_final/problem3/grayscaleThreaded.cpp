#include <iostream>
#include <vector>
#include <string>
#include <thread>
#include <chrono>
#include <cstring>
#include <png.h>

//threaded solution - dummy code
void grayscaleThreaded(png_bytep* image, int width, int height, int channels, int numThreads) {
    auto processRows = [&](int startY, int endY) {
        for (int y = startY; y < endY; ++y) {
            png_bytep row = image[y];
            for (int x = 0; x < width; ++x) {
                png_bytep px = &(row[x * channels]);
                uint8_t gray = static_cast<uint8_t>(0.3 * px[0] + 0.59 * px[1] + 0.11 * px[2]);
                px[0] = px[1] = px[2] = gray;
            }
        }
    };

    std::vector<std::thread> threads;
    int rowsPerThread = height / numThreads;
    int extraRows = height % numThreads;
    int startRow = 0;

    for (int i = 0; i < numThreads; ++i) {
        int endRow = startRow + rowsPerThread;
        if (i < extraRows) {
            endRow++;
        }
        threads.emplace_back(processRows, startRow, endRow);
        startRow = endRow;
    }

    for (auto& t : threads) {
        t.join();
    }
}
