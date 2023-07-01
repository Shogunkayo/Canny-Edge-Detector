#include "Image.h"

#include <cmath>
#include <chrono>

using namespace std;

int main(int argc, char** argv) {

	Image img("images/man.jpg");
    int img_size = img.w * img.h;

	// Convert to Grayscale
	img.grayscale_avg();
    img.write("test.jpg");

    cout << img.size << endl;

	Image gray_img(img.w, img.h, 1);
	for(uint64_t k=0; k<img_size; k++) {
		gray_img.data[k] = img.data[img.channels*k];
	}

    cout << gray_img.size << endl;
	gray_img.write("output/gray.png");

	// Apply Gaussian Blur
	Image blur_img(img.w, img.h, 1);
	double gauss[9] = {
		1/16.0, 2/16.0, 1/16.0,
		2/16.0, 4/16.0, 2/16.0,
		1/16.0, 2/16.0, 1/16.0
	};

	img.convolve_zero_pad(0, 3, 3, gauss, 1, 1);
	for(uint64_t k=0; k<img_size; ++k) {
		blur_img.data[k] = img.data[k];
	}
	blur_img.write("output/blur.png");


	// Edge detection
	double* tx = new double[img_size];
	double* ty = new double[img_size];
	double* gx = new double[img_size];
	double* gy = new double[img_size];

	for(uint32_t c=1; c < blur_img.w - 1; c++) {
		for(uint32_t r=0; r < blur_img.h; r++) {
			tx[r*blur_img.w + c] = blur_img.data[r*blur_img.w + c + 1] - blur_img.data[r*blur_img.w + c - 1];
			ty[r*blur_img.w + c] = 47*blur_img.data[r*blur_img.w + c + 1] + 162*blur_img.data[r*blur_img.w + c] + 47*blur_img.data[r*blur_img.w + c - 1];
		}
	}
	for(uint32_t c=1; c < blur_img.w - 1; c++) {
		for(uint32_t r=1; r < blur_img.h - 1; r++) {
			gx[r*blur_img.w + c] = 47*tx[(r+1)*blur_img.w + c] + 162*tx[r*blur_img.w + c] + 47*tx[(r-1)*blur_img.w + c];
			gy[r*blur_img.w + c] = ty[(r+1)*blur_img.w + c] - ty[(r-1)*blur_img.w + c];
		}
	}

	double max_x = -INFINITY, max_y = -INFINITY, min_x = INFINITY, min_y = INFINITY;
	for(uint64_t k=0; k<img_size; k++) {
		max_x = fmax(max_x, gx[k]);
		max_y = fmax(max_y, gy[k]);
		min_x = fmin(min_x, gx[k]);
		min_y = fmin(min_y, gy[k]);
	}
	Image Gx(img.w, img.h, 1);
	Image Gy(img.w, img.h, 1);
	for(uint64_t k=0; k < img_size; k++) {
		Gx.data[k] = (uint8_t)(255*(gx[k] - min_x) / (max_x - min_x));
		Gy.data[k] = (uint8_t)(255*(gy[k] - min_y) / (max_y - min_y));
	}
	Gx.write("output/Gx.png");
	Gy.write("output/Gy.png");

	double threshold = 0.01;
	double* g = new double[img_size];
	double* theta = new double[img_size];
	double x, y;

	for(uint64_t k=0; k < img_size; k++) {
		x = gx[k];
		y = gy[k];
		g[k] = sqrt(x*x + y*y);
		theta[k] = atan2(y, x);
	}

	double max = -INFINITY, min = INFINITY;
	for(uint64_t k=0; k < img_size; k++) {
		max = fmax(max, g[k]);
		min = fmin(min, g[k]);
	}
	Image G(img.w, img.h, 1);
	Image GT(img.w, img.h, 3);

	double h, s, l;
	double v;
	for(uint64_t k=0; k < img_size; k++) {
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

		G.data[k] = (uint8_t)(255 * v);
	}
	G.write("output/G.png");
	GT.write("output/GT.png");

    Image non_max(img.w, img.h, 1);
    for (uint64_t i=0; i < img_size; i++) {
        non_max.data[i] = 0;
    }

    for (uint64_t i=0; i < img_size; i++) {
        uint8_t q, r;
        
        // angle 0
        if ((0 <= theta[i] < 22.5) || (157.5 <= theta[i] <= 180)) {
            q = i + 1 < img_size ? G.data[i + 1] : 255;
            r = i - 1 > 0 ? G.data[i - 1] : 255;
        }

        // angle 45
        else if (22.5 <= theta[i] < 67.5) {
            q = i + img.w - 1 < img_size ? G.data[i + img.w - 1] : 255;
            r = i - img.w + 1 > 0 ? G.data[i - img.w + 1] : 255;
        }

        // angle 90
        else if (67.5 <= theta[i] < 112.5) {
            q = i + img.w < img_size ? G.data[i + img.w] : 255;
            r = i - img.w > 0 ? G.data[i - img.w] : 255;
        }

        // angle 135
        else if (112.5 <= theta[i] < 157.5) {
            q = i - img.w - 1 > 0 ? G.data[i - img.w - 1] : 255;
            r = i + img.w + 1 < img_size ? G.data[i + img.w + 1] : 255;
        }

        if ((G.data[i] >= q) && (G.data[i] >= r)) {
            non_max.data[i] = G.data[i];
        }
        else {
            non_max.data[i] = 0;
        }
    }

    non_max.write("output/non_max.png");

    // Double Threshold
    Image double_thresh (img.w, img.h, 1);

    double low_thresh_r = 0.05;
    double high_thresh_r = 0.09;
    
    max = -INFINITY;
    for (uint64_t i=0; i < img_size; i++) {
        max = fmax(max, non_max.data[i]);
    }

    double high_thresh = max * high_thresh_r;
    double low_thresh = high_thresh * low_thresh_r;

    uint8_t weak = 25;
    uint8_t strong = 255;

    for (uint64_t i=0; i < img_size; i++) {
        if (non_max.data[i] >= high_thresh) {
            double_thresh.data[i] = strong;
        }
        else if (non_max.data[i] >= low_thresh) {
            double_thresh.data[i] = weak;
        }
        else{
            double_thresh.data[i] = 0;
        }
    }

    double_thresh.write("output/double_thresh.png");

    // Edge Hysteresis
    Image canny(img.w, img.h, 1);
    for (uint64_t i=0; i < img_size; i++) {
        canny.data[i] = double_thresh.data[i];
    }

    for (uint64_t i=0; i < img_size; i++) {
        if (double_thresh.data[i] == weak) {
            if (
                    (i - img.w - 1 > 0 && double_thresh.data[i - img.w - 1] == strong) ||
                    (i - img.w > 0 && double_thresh.data[i - img.w] == strong) ||
                    (i - img.w + 1 > 0 && double_thresh.data[i - img.w + 1] == strong) ||
                    (i - 1 > 0 && double_thresh.data[i - 1] == strong) || 
                    (i + 1 > 0 && double_thresh.data[i + 1] == strong) ||
                    (i + img.w - 1 > 0 && double_thresh.data[i + img.w - 1] == strong) ||
                    (i + img.w > 0 && double_thresh.data[i + img.w] == strong) ||
                    (i + img.w + 1 > 0 && double_thresh.data[i + img.w + 1] == strong)
                ) 
                canny.data[i] = strong;
            else
                canny.data[i] = 0;
        }
    }

    canny.write("output/canny.png");

    delete[] tx;
    delete[] ty;
	delete[] gx;
	delete[] gy;
	delete[] g;
	delete[] theta;

	return 0;
}
