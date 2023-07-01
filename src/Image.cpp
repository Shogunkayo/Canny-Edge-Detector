#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#define BYTE_BOUND(value) value < 0 ? 0 : (value > 255 ? 255 : value)

#include "Image.h"
#include "stb_image.h"
#include "stb_image_write.h"


Image::Image(const char* filename) {
    if(read(filename)) {
        std::cout << "Read file: " << filename << std::endl;
        size = w * h * channels;
    }

    else {
        std::cout << "Read failed" << std::endl;
    }
}


Image::Image(int w, int h, int channels) : w(w), h(h), channels(channels) {
    //We are initializing the data members of an object from Image.h
    
    size = w * h * channels;
    data = new uint8_t[size];
}


Image::Image(const Image& img) : Image(img.w, img.h, img.channels) {
    memcpy(data, img.data, size);
}


Image::~Image() {
    stbi_image_free(data);  // Just cleans up all the data
}


bool Image::read(const char* filename) {
    // This function is used in the constructor
    
    // (Filename, Width, Height, Number of channels, no of channels to force)
    data = stbi_load(filename, &w, &h, &channels, 0);
    return data != NULL;
}


bool Image::write(const char* filename) {
    ImageType type = getFileType(filename);

    //Stores different codes for success of writing different file types
    int success;

    switch(type) {
        case PNG:
            success = stbi_write_png(filename, w, h, channels, data, w * channels);
            break;
        
        case BMP:
            success = stbi_write_bmp(filename, w, h, channels, data);
            break;

        case JPG:
            success = stbi_write_jpg(filename, w, h, channels, data, 100);
            break;

        case TGA:
            success = stbi_write_tga(filename, w, h, channels, data);
            break;
    }

    return success != 0;
}


ImageType Image::getFileType(const char* filename) {
    //Gives the extension of the file
    const char* ext = strrchr(filename, '.');

    if(ext != nullptr) {
        if(strcmp(ext, ".png") == 0) {
            return PNG;
        }

        else if(strcmp(ext, ".jpg") == 0) {
            return JPG;
        }

        else if(strcmp(ext, ".bmp") == 0) {
            return BMP;
        }

        else if(strcmp(ext, ".tga") == 0) {
            return TGA;
        }

        else {
            return PNG; //Default case
        }
    }
}

Image& Image::grayscale_avg() {
    if (channels < 3) {
        std::cout << "Expected 3 channels, found " << channels << std::endl;
    }
    else {
        for (int i=0; i < size; i += channels) {
            int gray = (data[i] + data[i+1] + data[i+2])/3;
            memset(data+i, gray, 3);
        }

        uint8_t *temp = new uint8_t[h*w];
        for (uint64_t i=0; i < size; i++) {
            temp[i] = data[channels * i];
        }

        data = new uint8_t[h*w];
        for (uint64_t i=0; i < h*w; i++) {
            data[i] = temp[i];
        }

        delete[] temp;
    }

    return *this;
}

Image& Image::grayscale_lum() {
    // Weighted average

    if (channels < 3) {
        std::cout << "Expected 3 channels, found " << channels << std::endl;
    }
    else {
        for (int i=0; i < size; i += channels) {
            int gray = 0.2126*data[i] + 0.7152*data[i+1] + 0.0722*data[i+2];
            memset(data+i, gray, 3);
        }
    }

    return *this;
}

Image& Image::colorMask(float r, float g, float b) {
    if (channels < 3) {
        std::cout << "Expected 3 channels, found " << channels << std::endl;
    }

    else {
        for (int i=0; i < size; i += channels) {
            data[i] *= r;
            data[i+1] *= g;
            data[i+2] *= b;
        }
    }

    return *this;
}

Image& Image::convolve_border(uint8_t channel, uint32_t ker_w, 
        uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc) {

    uint8_t sum[w*h];
    uint64_t center = cr * ker_w + cc;

    for (uint64_t k=channel; k < size; k += channels) {
        double c = 0;
        for (long i = -(long)cr; i < (long)ker_h - cr; i++) {
            long row = ((long)k / channels) / w-i;

            if (row < 0) {
                row = 0;
            }
            else if (row > h-1) {
                row = h - 1;
            }

            for (long j = -(long)cc; j < (long)ker_w - cc; j++) {
                long col = ((long)k / channels) % w-j;

                if (col < 0) {
                    col = 0;
                }
                else if (col > w-1) {
                    col = w - 1;
                }
                c += ker[center + i * (long)ker_w + j] 
                    * data[(row*w+col)*channels + channel];
            }
        }
        
        sum[k/channels] = (uint8_t)BYTE_BOUND(round(c));
    }

    for (uint64_t k=channel; k < size; k += channels) {
        data[k] = sum[k/channels];
    }

    return *this;
}

Image& Image::convolve_zero_pad(uint8_t channel, uint32_t ker_w, 
        uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc) {

    uint8_t sum[w*h];
    uint64_t center = cr * ker_w + cc;

    for (uint64_t k=channel; k < size; k += channels) {
        double c = 0;
        for (long i = -(long)cr; i < (long)ker_h - cr; i++) {
            long row = ((long)k / channels) / w-i;

            if (row < 0 || row > h-1) {
                continue;
            }

            for (long j = -(long)cc; j < (long)ker_w - cc; j++) {
                long col = ((long)k / channels) % w-j;

                if (col < 0 || col > w-1) {
                    continue;
                }
                c += ker[center + i * (long)ker_w + j] 
                    * data[(row*w+col)*channels + channel];
            }
        }
        
        sum[k/channels] = (uint8_t)BYTE_BOUND(round(c));
    }

    for (uint64_t k=channel; k < size; k += channels) {
        data[k] = sum[k/channels];
    }

    return *this;
}

Image& Image::convolve_cyclic(uint8_t channel, uint32_t ker_w, 
        uint32_t ker_h, double ker[], uint32_t cr, uint32_t cc) {

    uint8_t sum[w*h];
    uint64_t center = cr * ker_w + cc;

    for (uint64_t k=channel; k < size; k += channels) {
        double c = 0;
        for (long i = -(long)cr; i < (long)ker_h - cr; i++) {
            long row = ((long)k / channels) / w-i;

            if (row < 0) {
                row = row % h + h;
            }
            else if (row > h-1) {
                row %= h;
            }

            for (long j = -(long)cc; j < (long)ker_w - cc; j++) {
                long col = ((long)k / channels) % w-j;

                if (col < 0) {
                    col = col % w + w;
                }
                else if (col > w-1) {
                    col %= w;
                }
                c += ker[center + i * (long)ker_w + j] 
                    * data[(row*w+col)*channels + channel];
            }
        }
        
        sum[k/channels] = (uint8_t)BYTE_BOUND(round(c));
    }

    for (uint64_t k=channel; k < size; k += channels) {
        data[k] = sum[k/channels];
    }

    return *this;
}

Image& Image::flipY() {
    uint8_t temp[4];
    uint8_t* px1;
    uint8_t* px2;

    for (int y=0; y < h; y++) {
        for (int x = 0; x < w/2; x++) {
            px1 = &data[(x + y * w) * channels];
            px2 = &data[((w - 1 - x) + y * w) * channels];

            memcpy(temp, px1, channels);
            memcpy(px1, px2, channels);
            memcpy(px2, temp, channels);
        }
    }

    return *this;
}

Image& Image::flipX() {
    uint8_t temp[4];
    uint8_t* px1;
    uint8_t* px2;

    for (int x=0; x < w; x++) {
        for (int y=0; y < h/2; y++) {
            px1 = &data[(x + y * w) * channels];
            px2 = &data[(x + (h - 1 - y) * w) * channels];

            memcpy(temp, px1, channels);
            memcpy(px1, px2, channels);
            memcpy(px2, temp, channels);
        }
    }

    return *this;
}

Image& Image::crop(uint16_t cx, uint16_t cy, uint16_t cw, uint16_t ch) {
    uint8_t* cropped = new uint8_t[cw * ch * channels];
    memset(cropped, 0, cw * ch * channels);

    for (uint16_t y = 0; y < ch; y++) {
        if (y + cy >= h) {
            break;
        }

        for (uint16_t x = 0; x < cw; x++) {
            if (x + cx >= w) {
                break;
            }
            memcpy(&cropped[(x + y * cw) * channels], 
                    &data[(x + cx + (y + cy) * w) * channels ], channels);
        }
    }

    w = cw;
    h = ch;
    size = w * h * channels;

    delete[] data;
    data = cropped;
    cropped = nullptr;

    return *this;
}

Image& Image::blur_gaussian() {

    

    return *this;
}

Image& Image::scharr(double threshold=0.09) {
    return *this;
}

Image& Image::sobel(double threshold=0.09) {
    return *this;
}

Image& Image::non_max() {
    return *this;
}

Image& Image::double_threshold(double low_r=0.05, double high_r=0.09, 
        uint8_t strong=255, uint8_t weak=25) {
    return *this;
}

Image& Image::edge_hysteresis(uint8_t strong=255, uint8_t weak=25) {
    return *this;
}

Image& Image::canny() {
    return *this;
}
