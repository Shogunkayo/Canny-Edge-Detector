#include "Image.h"

#include <cstdlib>
#include <cmath>
#include <chrono>


int main(int argc, char** argv) {

	Image img("images/man.jpg");
    img.canny();

	return 0;
}
