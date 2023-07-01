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
        for (uint64_t k=0; k < h*w; k++) {
            temp[k] = data[k * channels];
        }

        delete[] data;
        data = new uint8_t[h*w];
        for (uint64_t i=0; i < h*w; i++) {
            data[i] = temp[i]; 
        }

        delete[] temp;
        channels = 1;
        size = h*w;

        write("output/grayscale_avg.png");
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
    double gauss[9] = {
		1/16.0, 2/16.0, 1/16.0,
		2/16.0, 4/16.0, 2/16.0,
		1/16.0, 2/16.0, 1/16.0
	};

    convolve_zero_pad(0, 3, 3, gauss, 1, 1);
    write("output/blur_gaussian.png");

    return *this;
}

Image& Image::scharr(double threshold) {
    grayscale_avg();
    blur_gaussian();

    double* tx = new double[size];
	double* ty = new double[size];
	double* gx = new double[size];
	double* gy = new double[size];

	for(uint32_t c=1; c < w - 1; c++) {
		for(uint32_t r=0; r < h; r++) {
			tx[r*w + c] = data[r*w + c + 1] - data[r*w + c - 1];
			ty[r*w + c] = 47*data[r*w + c + 1] + 162*data[r*w + c] + 47*data[r*w + c - 1];
		}
	}
	for(uint32_t c=1; c < w - 1; c++) {
		for(uint32_t r=1; r < h - 1; r++) {
			gx[r*w + c] = 47*tx[(r+1)*w + c] + 162*tx[r*w + c] + 47*tx[(r-1)*w + c];
			gy[r*w + c] = ty[(r+1)*w + c] - ty[(r-1)*w + c];
		}
	}

	double max_x = -INFINITY, max_y = -INFINITY, min_x = INFINITY, min_y = INFINITY;
	for(uint64_t k=0; k < size; k++) {
		max_x = fmax(max_x, gx[k]);
		max_y = fmax(max_y, gy[k]);
		min_x = fmin(min_x, gx[k]);
		min_y = fmin(min_y, gy[k]);
	}

	Image Gx(w, h, 1);
	Image Gy(w, h, 1);
	for(uint64_t k=0; k < size; k++) {
		Gx.data[k] = (uint8_t)(255*(gx[k] - min_x) / (max_x - min_x));
		Gy.data[k] = (uint8_t)(255*(gy[k] - min_y) / (max_y - min_y));
	}
	Gx.write("output/Gx.png");
	Gy.write("output/Gy.png");

	double* g = new double[size];
	theta = new double[size];
	double x, y;

	for(uint64_t k=0; k < size; k++) {
		x = gx[k];
		y = gy[k];
		g[k] = sqrt(x*x + y*y);
		theta[k] = atan2(y, x);
	}

	double max = -INFINITY, min = INFINITY;
	for(uint64_t k=0; k < size; k++) {
		max = fmax(max, g[k]);
		min = fmin(min, g[k]);
	}
	Image G(w, h, 1);
	Image GT(w, h, 3);

	double h, s, l;
	double v;
	for(uint64_t k=0; k < size; k++) {
		// Hue = theta between 0 and 360
		h = theta[k] * 180.0/M_PI + 180.;

		//v is the relative edge strength
		if(max == min) {
			v = 0;
		}
		else {
			v = (g[k] - min) / (max - min);
		}
		s = l = v;

		// HSL to RGB
		double c = (1 - abs(2*l - 1)) * s;
		double x = c * (1 - abs(fmod((h/60), 2) - 1));
		double m = l - c/2;

		double r = 0, g = 0, b = 0;
		if(h < 60) {
			r = c;
			g = x;
		}
		else if(h < 120) {
			r = x;
			g = c;
		}
		else if(h < 180) {
			g = c;
			b = x;
		}
		else if(h < 240) {
			g = x;
			b = c;
		}
		else if(h < 300) {
			b = c;
			r = x;
		}
		else {
			b = x;
			r = c;
		}

		uint8_t red, green, blue;
		red = (uint8_t)(255 *(r + m));
		green = (uint8_t)(255 *(g + m));
		blue = (uint8_t)(255 *(b + m));

		GT.data[k*3] = red;
		GT.data[k*3 + 1] = green;
		GT.data[k*3 + 2] = blue;

		data[k] = (uint8_t)(255 * v);
	}
	write("output/G.png");
	GT.write("output/GT.png");

    delete[] tx;
    delete[] ty;
    delete[] gx;
    delete[] gy;
    delete[] g;

    return *this;
}

Image& Image::sobel(double threshold) {
    return *this;
}

Image& Image::non_max() {
    if (theta == NULL) {
        std::cout << "Requires theta" << std::endl;
    }
    else {
        for (uint64_t i=0; i < size; i++) {
            uint8_t q, r;

            // angle 0
            if ((0 <= theta[i] < 22.5) || (157.5 <= theta[i] <= 180)) {
                q = i + 1 < size ? data[i + 1] : 255;
                r = i - 1 > 0 ? data[i - 1] : 255;
            }

            // angle 45
            else if (22.5 <= theta[i] < 67.5) {
                q = i + w - 1 < size ? data[i + w - 1] : 255;
                r = i - w + 1 > 0 ? data[i - w + 1] : 255;
            }

            // angle 90
            else if (67.5 <= theta[i] < 112.5) {
                q = i + w < size ? data[i + w] : 255;
                r = i - w > 0 ? data[i - w] : 255;
            }

            // angle 135
            else if (112.5 <= theta[i] < 157.5) {
                q = i - w - 1 > 0 ? data[i - w - 1] : 255;
                r = i + w + 1 < size ? data[i + w + 1] : 255;
            }

            if ((data[i] >= q) && (data[i] >= r)) {
                data[i] = data[i];
            }
            else {
                data[i] = 0;
            }
        }

        write("output/non_max.png");

    }


    return *this;
}

Image& Image::double_threshold(double low_r, double high_r, uint8_t strong, 
        uint8_t weak) {

    double max = -INFINITY;
    for (uint64_t i=0; i < size; i++) {
        max = fmax(max, data[i]);
    }

    double high_thresh = max * high_r;
    double low_thresh = high_thresh * low_r;

    for (uint64_t i=0; i < size; i++) {
        if (data[i] >= high_thresh) {
            data[i] = strong;
        }
        else if (data[i] >= low_thresh) {
            data[i] = weak;
        }
        else{
            data[i] = 0;
        }
    }

    write("output/double_thresh.png");

    return *this;
}

Image& Image::edge_hysteresis(uint8_t strong, uint8_t weak) {
    uint8_t* temp = new uint8_t[size];
    for (uint64_t i=0; i < size; i++) {
        temp[i] = data[i];
    }

    for (uint64_t i=0; i < size; i++) {
        if (data[i] == weak) {
            if (
                    (i - w - 1 > 0 && data[i - w - 1] == strong) ||
                    (i - w > 0 && data[i - w] == strong) ||
                    (i - w + 1 > 0 && data[i - w + 1] == strong) ||
                    (i - 1 > 0 && data[i - 1] == strong) || 
                    (i + 1 > 0 && data[i + 1] == strong) ||
                    (i + w - 1 > 0 && data[i + w - 1] == strong) ||
                    (i + w > 0 && data[i + w] == strong) ||
                    (i + w + 1 > 0 && data[i + w + 1] == strong)
                ) 
                temp[i] = strong;
            else
                temp[i] = 0;
        }
    }

    for (uint64_t i=0; i < size; i++) {
        data[i] = temp[i];
    }

    delete[] temp;
    return *this;
}

Image& Image::canny() {
    scharr();
    non_max();
    double_threshold();
    edge_hysteresis();

    write("output/canny.png");

    return *this;
}
