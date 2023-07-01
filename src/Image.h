#include <stdint.h> //Contains uint8_t and others
#include <iostream>

enum ImageType {
    //An enum to store and access our various image file types

    PNG, JPG, BMP, TGA
};

struct Image {
    // Data members
    uint8_t * data = NULL;   // Contains all picture data in 1 byte
    size_t size = 0;    // Keeps track of size of data store
    int w;  // Width of image
    int h;  // Height of image
    int channels;   // Number of channels of image
    double* theta = NULL;

    // Constructors and Destructors
    Image(const char*);  //a Constructor that takes in only the file name
    Image(int w, int h, int channels);  //Constructor that takes in h, w and channels
    Image(const Image& img); // Copy constructor
    ~Image();   // Destructor

    // Functions
    bool read(const char* filename);    // Reads in from a file
    bool write(const char* filename);   // Writes into the file 
    ImageType getFileType(const char* filename);    // Returns image file type

    // -- Grayscale
    Image& grayscale_avg();
    Image& grayscale_lum();

    // -- Colormask
    Image& colorMask(float r, float g, float b);

    // -- Flipping
    Image& flipX();
    Image& flipY();

    // -- Crop
    Image& crop(uint16_t cx, uint16_t cy, uint16_t cw, uint16_t ch);

    // -- Convolution
    Image& convolve_zero_pad(uint8_t channel, uint32_t ker_w, uint32_t ker_h,
            double ker[], uint32_t cr, uint32_t cc);
    Image& convolve_border(uint8_t channel, uint32_t ker_w, uint32_t ker_h,
            double ker[], uint32_t cr, uint32_t cc);
    Image& convolve_cyclic(uint8_t channel, uint32_t ker_w, uint32_t ker_h,
            double ker[], uint32_t cr, uint32_t cc);

    // -- Edges
    Image& blur_gaussian();
    Image& scharr(double threshold=0.09);
    Image& sobel(double threshold=0.09);
    Image& non_max();
    Image& double_threshold(double low_r=0.05, double high_r=0.09, uint8_t strong=255, uint8_t weak=25);
    Image& edge_hysteresis(uint8_t strong=255, uint8_t weak=25);
    Image& canny();
};
