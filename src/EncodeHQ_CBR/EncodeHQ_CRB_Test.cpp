/*********************************************************************/
/* EncodeHQCBR.cpp                                                   */
/* Author: Tim Borer                                                 */
/* This version 18th September 2013                                  */
/*                                                                   */
/* Reads image data in from a planar file.                           */
/* Compresses image using VC-2 High Quality profile @CBR.            */
/* Write compressed transform data out (not complete stream).        */
/* It is not necessarily complet nor korrect.                        */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file   */
/*********************************************************************/

const char version[] = __DATE__ " @ " __TIME__;
const char summary[] = "Encodes an uncompressed planar video file with VC-2 High Quality profile at constant bit rate";
const char description[] = "\
This program compresses an image sequence using SMPTE VC-2 HQ profile.\n\
It implements constant bit rate coding.\n\
The bit rate is specified by defining the number of compressed bytes per frame.\n\
Its primary output is the compressed bytes. However it may produce alternative outputs which are:\n\
  1 the wavelet transform of the input\n\
  2 the quantised wavelet coefficients\n\
  3 the quantisation indices used for each slice\n\
  4 compressed bytes\n\
  5 VC2 bitstream (default output)\n\
  6 the decoded sequence\n\
  7 the PSNR for each frame\n\
Input and output (where appropriate) are in planar format (4:4:4, 4:2:2, 4:2:0 or RGB).\n\
There can be 1 to 4 bytes per sample and the data is left (MSB) justified.\n\
Data is assumed offset binary (which is fine for both YCbCr or RGB).\n\
\n\
Example: EncodeHQ-CBR -v -x 1920 -y 1080 -f 4:2:2 -l 10 -k LeGall -d 3 -u 1 -a 2 -s 829440 -i inFileName outFileName";
const char* details[] = { version, summary, description };

#include <cstdlib> //for EXIT_SUCCESS, EXIT_FAILURE, atoi
#include <stdexcept> //For standard logic errors
#include <iostream> //For cin, cout, cerr
#include <string>
#include <fstream>
#include <cstdio> // for perror
#include <iomanip> // For reporting stats only
#include <algorithm>
#include <functional>
#include <cmath>

#include "EncodeParams.h"
#include "Arrays.h"
#include "Picture.h"
#include "Frame.h"
#include "WaveletTransform.h"
#include "Quantisation.h"
#include "Slices.h"
#include "DataUnit.h"
#include "Utils.h"

using std::cout;
using std::cin;
using std::cerr;
using std::clog;
using std::endl;
using std::string;
using std::ifstream;
using std::ofstream;
using std::ios_base;
using std::filebuf;
using std::streambuf;
using std::fill_n;

using arrayio::ioFormat;  // enu
using arrayio::format;    // class
using arrayio::wordWidth;  //class
using arrayio::bitDepth;   //class
using arrayio::offset;    //class
using arrayio::left_justified;
using arrayio::right_justified;

// Calculate quantisation indices using a binary search
const Array2D quantIndices(const Picture& coefficients,
	const Array1D& qMatrix,
	const Array2D& sliceBytes,
	const int scalar) {
	const int ySlices = sliceBytes.shape()[0];
	const int xSlices = sliceBytes.shape()[1];
	// Create an empty array of indices to fill and return
	Array2D indices(extents[ySlices][xSlices]);
	// Wavelet depth & number of subbands derived from dimensions of qMatrix
	const int numberOfSubbands = qMatrix.size();
	const int waveletDepth = (numberOfSubbands - 1) / 3;
	const PictureArray slices = split_into_blocks(coefficients, ySlices, xSlices);
	for (int row = 0; row<ySlices; ++row) {
		for (int column = 0; column<xSlices; ++column) {
			// Available bytes is the size of slice less 4 byte overhead
			const int bytesAvailable = sliceBytes[row][column] - 4;
			int trialQ = 63;
			int q = 127;
			int delta = 64;
			while (delta>0) {
				delta >>= 1;
				const Picture trialSlice = quantise_transform_np(slices[row][column], trialQ, qMatrix);
				int bytesRequired = component_slice_bytes(trialSlice.y(), waveletDepth, scalar);
				bytesRequired += component_slice_bytes(trialSlice.c1(), waveletDepth, scalar);
				bytesRequired += component_slice_bytes(trialSlice.c2(), waveletDepth, scalar);
				if (bytesRequired <= bytesAvailable) {
					if (trialQ<q) q = trialQ;
					trialQ -= delta;
				}
				else {
					trialQ += delta;
				}
			}
			indices[row][column] = q;
		}
	}
	return indices;
}

int main(void) {



		  //Open input file in binary mode.
		int bits = 8;
		ColourFormat chromaFormat = CF422;  // {UNKNOWN, CF444, CF422, CF420, RGB};

		string inFileName = "D:/resource/history/TEST8.PPM";
		string outFileName = "Haar1D1YUV422Bits10.drc";

		const WaveletKernel kernel = LeGall; // {DD97, LeGall, DD137, Haar0, Haar1, Fidelity, Daub97, NullKernel};
		const int waveletDepth = 3;
		float CompressedRate = 2.0;   // 2, 5, 8
		int height;
		int width;
		int MaxValue;
		const int sliceScalar = 1;

		int nbytes;
		if (bits <= 8)
		{
			nbytes = 1;
		}
		else
		{
			nbytes = 2;
		}

		int ySize;
		int xSize;

		ifstream input;
		input.open(inFileName.c_str(), ios_base::in | ios_base::binary);
		if (!input)
		{
			cerr << "Error: failed to open input file " << inFileName << endl;
			return 0;
		}

		input >> wordWidth(nbytes);
		input >> bitDepth(bits);
		input >> right_justified;
		input >> offset(0);

		string headP6;
		input >> headP6;
		input >> height;
		input >> width;
		input >> MaxValue;
		input.get();

		int compressedBytes;
		switch (chromaFormat) {
		case RGB:
		case CF444:
			compressedBytes = floor(height*width*nbytes*3.0 / CompressedRate);
			ySize = 1;
			xSize = 1;
			break;
		case CF422:
			compressedBytes = floor(height*width*nbytes*2.0 / CompressedRate);
			ySize = 1;
			xSize = 2;
			break;
		case CF420:
			compressedBytes = floor(height*width*nbytes*1.5 / CompressedRate);
			ySize = 2;
			xSize = 2;
			break;
		}
		//Read pixel data, line by line (to maximise cache occupancy).

		const Shape2D  ppmsize = { { height, 3 * width } };
		Array2D  RGBArray(ppmsize);

		input >> RGBArray;
		PictureFormat pctFormat(height, width, chromaFormat);

		const Shape2D  result = { { height, width } };
		const Shape2D  UVresult422 = { { height, width / 2 } };
		const Shape2D  UVresult420 = { { height / 2, width / 2 } };

		Array2D  YArray;
		Array2D  UArray;
		Array2D  VArray;


		const int UVHeight = height + 2;
		const int UVWidth = width + 2;
		const int UVImageSize = UVHeight*UVWidth;

		int *ULine = (new int[width + 2]) + 1;
		int *VLine = (new int[width + 2]) + 1;
		int *UImage = (new int[UVImageSize]) + UVWidth + 1;
		int *VImage = (new int[UVImageSize]) + UVWidth + 1;
		ULine[-1] = ULine[width] = 128;
		VLine[-1] = VLine[width] = 128;


		int Y, U, V;
		int R, G, B;
		switch (chromaFormat) {
		case RGB:
			YArray.resize(result);
			UArray.resize(result);
			VArray.resize(result);
			for (int line = 0; line<height; line++) {
				for (register int pixel = 0; pixel<width; ++pixel) {
					R = RGBArray[line][3 * pixel];
					G = RGBArray[line][3 * pixel + 1];
					B = RGBArray[line][3 * pixel + 2];

					YArray[line][pixel] = R;
					UArray[line][pixel] = G;
					VArray[line][pixel] = B;
				}
			}
			break;
		case CF444:
			YArray.resize(result);
			UArray.resize(result);
			VArray.resize(result);
			for (int line = 0; line<height; line++) {
				for (register int pixel = 0; pixel<width; ++pixel) {
					R = RGBArray[line][3 * pixel];
					G = RGBArray[line][3 * pixel + 1];
					B = RGBArray[line][3 * pixel + 2];

					//Convert RGB to YUV
					Y = ((66 * R + 129 * G + 25 * B + 128) >> 8) + 16;
					U = ((-38 * R - 74 * G + 112 * B + 128) >> 8) + 128;
					V = ((112 * R - 94 * G - 18 * B + 128) >> 8) + 128;

					//Clip YUV values
					YArray[line][pixel] = static_cast<int>((Y<0) ? 0 : ((Y>MaxValue) ? MaxValue : Y));
					UArray[line][pixel] = static_cast<int>((U<0) ? 0 : ((U>MaxValue) ? MaxValue : U));
					VArray[line][pixel] = static_cast<int>((V<0) ? 0 : ((V>MaxValue) ? MaxValue : V));
				}
			}
			break;
		case CF422:

			YArray.resize(result);
			UArray.resize(UVresult422);
			VArray.resize(UVresult422);

			for (int line = 0; line<height; ++line) {
				for (int pixel = 0; pixel<width; ++pixel) {

					R = RGBArray[line][3 * pixel];
					G = RGBArray[line][3 * pixel + 1];
					B = RGBArray[line][3 * pixel + 2];

					//Convert RGB to YUV
					Y = ((66 * R + 129 * G + 25 * B + 128) >> 8) + 16;
					U = ((-38 * R - 74 * G + 112 * B + 128) >> 8) + 128;
					V = ((112 * R - 94 * G - 18 * B + 128) >> 8) + 128;

					//Clip Y ready for output & copy UV ready for filtering
					YArray[line][pixel] = static_cast<int>((Y<0) ? 0 : ((Y>MaxValue) ? MaxValue : Y));
					ULine[pixel] = U;
					VLine[pixel] = V;
				}

				for (int pixel = 0; pixel<width; pixel += 2) {

					//Filter line
					U = ((ULine[pixel - 1] + 2 * ULine[pixel] + ULine[pixel + 1] + 2) >> 2);
					V = ((VLine[pixel - 1] + 2 * VLine[pixel] + VLine[pixel + 1] + 2) >> 2);

					//Clip and copy UV to output buffer
					UArray[line][pixel / 2] = static_cast<int>((U<0) ? 0 : ((U>MaxValue) ? MaxValue : U));
					VArray[line][pixel / 2] = static_cast<int>((V<0) ? 0 : ((V>MaxValue) ? MaxValue : V));
				}
			}
			break;
		case CF420:
			YArray.resize(result);
			UArray.resize(UVresult420);
			VArray.resize(UVresult420);

			fill_n(&UImage[-(UVWidth + 1)], UVImageSize, 128);
			fill_n(&VImage[-(UVWidth + 1)], UVImageSize, 128);
			for (int line = 0; line<height; ++line) {
				for (int pixel = 0; pixel<width; ++pixel) {

					R = RGBArray[line][3 * pixel];
					G = RGBArray[line][3 * pixel + 1];
					B = RGBArray[line][3 * pixel + 2];

					//Convert RGB to YUV
					Y = ((66 * R + 129 * G + 25 * B + 128) >> 8) + 16;
					U = ((-38 * R - 74 * G + 112 * B + 128) >> 8) + 128;
					V = ((112 * R - 94 * G - 18 * B + 128) >> 8) + 128;

					//Clip Y ready for output & copy UV ready for filtering
					YArray[line][pixel] = static_cast<int>((Y<0) ? 0 : ((Y>MaxValue) ? MaxValue : Y));
					ULine[pixel] = U;
					VLine[pixel] = V;
				}

				for (int pixel = 0; pixel<width; pixel += 2) {
					//Filter line & horizontally subsample by 2
					//                UImage[line*UVWidth+pixel] = ULine[pixel];
					//                VImage[line*UVWidth+pixel] = VLine[pixel];
					UImage[line*UVWidth + pixel] = ((ULine[pixel - 1] + 2 * ULine[pixel] + ULine[pixel + 1] + 2) >> 2);
					VImage[line*UVWidth + pixel] = ((VLine[pixel - 1] + 2 * VLine[pixel] + VLine[pixel + 1] + 2) >> 2);
				}
			}

			for (int line = 0; line<height; line += 2) {
				for (int pixel = 0; pixel<width; pixel += 2) {

					//Filter vertically and subsample by 2
					//                U = UImage[line*UVWidth+pixel];
					//                V = VImage[line*UVWidth+pixel];
					U = ((UImage[(line - 1)*UVWidth + pixel] +
						2 * UImage[line*UVWidth + pixel] + UImage[(line + 1)*UVWidth + pixel] + 2) >> 2);
					V = ((VImage[(line - 1)*UVWidth + pixel] +
						2 * VImage[line*UVWidth + pixel] + VImage[(line + 1)*UVWidth + pixel] + 2) >> 2);

					//Clip and copy UV to output buffer
					UArray[line / 2][pixel / 2] = static_cast<int>((U<0) ? 0 : ((U>MaxValue) ? MaxValue : U));
					VArray[line / 2][pixel / 2] = static_cast<int>((V<0) ? 0 : ((V>MaxValue) ? MaxValue : V));
				}
			}
			break;
		default:
			cerr << "output Image format ERROR" << endl;
			return EXIT_FAILURE;
		} // switch


		delete[](&VLine[-1]);
		delete[](&ULine[-1]);
		delete[](&VImage[-(UVWidth + 1)]);
		delete[](&UImage[-(UVWidth + 1)]);


		Picture picture(pctFormat, YArray, UArray, VArray);

		const bool verbose = 1;
		const int lumaDepth = bits;
		int chromaDepth = bits;

		const bool interlaced = 0;
		const bool topFieldFirst = 0;
		const Output output = STREAM;

		
		ofstream outStream;
		outStream.open(outFileName.c_str(), ios_base::out | ios_base::binary);
		if (!outStream)
		{
			cerr << "Error: failed to open input file " << outFileName << endl;
			return 0;
		}

		PictureFormat format(height, width, chromaFormat);

		if (verbose) {
			clog << "bytes per sample= " << nbytes << endl;
			clog << "luma depth (bits) = " << lumaDepth << endl;
			clog << "chroma depth (bits) = " << chromaDepth << endl;
			clog << "height = " << format.lumaHeight() << endl;
			clog << "width = " << format.lumaWidth() << endl;
			clog << "chroma format = " << format.chromaFormat() << endl;
			clog << "interlaced = " << std::boolalpha << interlaced << endl;
			if (interlaced) clog << "top field first = " << std::boolalpha << topFieldFirst << endl;
			clog << "wavelet kernel = " << kernel << endl;
			clog << "wavelet depth = " << waveletDepth << endl;
			clog << "vertical slice size (in units of 2**(wavelet depth)) = " << ySize << endl;
			clog << "horizontal slice size (in units of 2**(wavelet depth)) = " << xSize << endl;
			clog << "compressed bytes = " << compressedBytes << endl;
			clog << "output = " << output << endl;
		}

		// Calculate number of slices per picture
		const int yTransformSize = ySize*utils::pow(2, waveletDepth);
		const int xTransformSize = xSize*utils::pow(2, waveletDepth);
		const int pictureHeight = ((interlaced) ? height / 2 : height);
		const int paddedPictureHeight = paddedSize(pictureHeight, waveletDepth);
		const int paddedWidth = paddedSize(width, waveletDepth);
		const int ySlices = paddedPictureHeight / yTransformSize;
		const int xSlices = paddedWidth / xTransformSize;
		if (paddedPictureHeight != (ySlices*yTransformSize)) {
			throw std::logic_error("Padded picture height is not divisible by slice height");
			return EXIT_FAILURE;
		}
		if (paddedWidth != (xSlices*xTransformSize)) {
			throw std::logic_error("Padded width is not divisible by slice width");
			return EXIT_FAILURE;
		}

		if (verbose) {
			clog << "Vertical slices per picture          = " << ySlices << endl;
			clog << "Horizontal slices per picture        = " << xSlices << endl;
			// Calculate slice bytes numerator and denominator
			const utils::Rational sliceBytesNandD =
				utils::rationalise((interlaced ? compressedBytes / 2 : compressedBytes), (ySlices*xSlices));
			const int SliceBytesNum = sliceBytesNandD.numerator;
			const int SliceBytesDenom = sliceBytesNandD.denominator;
			clog << "Slice bytes numerator                = " << SliceBytesNum << endl;
			clog << "Slice bytes denominator              = " << SliceBytesDenom << endl;
		}

		// Calculate the quantisation matrix
		const Array1D qMatrix = quantMatrix(kernel, waveletDepth);
		if (verbose) {
			clog << "Quantisation matrix = " << qMatrix[0];
			for (unsigned int i = 1; i<qMatrix.size(); ++i) {
				clog << ", " << qMatrix[i];
			}
			clog << endl;
		}


		//Forward wavelet transform
		if (verbose) clog << "Forward transform" << endl;
		Picture transform = waveletTransform(picture, kernel, waveletDepth);

		if (output == TRANSFORM) {
			//Write transform output as 4 byte 2's comp values
			clog << "Writing transform coefficients to output file" << endl;
			outStream << pictureio::wordWidth(4); //4 bytes per sample
			outStream << pictureio::signed_binary;
			outStream << transform;
			if (!outStream) {
				cerr << "Failed to write output file \"" << outFileName << "\"" << endl;
				return EXIT_FAILURE;
			}
		}

		// Choose quantisation indices to achieve a size of compressedBytes for the frame
		if (verbose) clog << "Determine quantisation indices" << endl;
		const int pictureBytes = (interlaced ? compressedBytes / 2 : compressedBytes);
		// Calculate number of bytes for each slice
		const Array2D bytes = slice_bytes(ySlices, xSlices, pictureBytes, sliceScalar);
		Array2D qIndices = quantIndices(transform, qMatrix, bytes, sliceScalar);

		if (verbose) clog << "Quantise transform coefficients" << endl;
		const Picture quantisedSlices = quantise_transform_np(transform, qIndices, qMatrix);


		// Split transform into slices
		if (verbose) clog << "Split quantised coefficients into slices" << endl;
		const PictureArray slices = split_into_blocks(quantisedSlices, ySlices, xSlices);

		int frame = 1;
		if (output == STREAM) {
			const int slicePrefix = 0;
			const Slices outSlices(slices, waveletDepth, qIndices);
			const WrappedPicture outWrapped(frame,
									kernel,
									waveletDepth,
									xSlices,
									ySlices,
									slicePrefix,
									sliceScalar,
									outSlices);

			//Write packaged output
			if (verbose) clog << "Writing compressed output to file" << endl;
			//for (int r = 0; r < bytes.shape()[0]; r++) {
			//	for (int c = 0; c < bytes.shape()[1]; c++)
			//		std::cout << bytes[r][c] << "\t";
			//	std::cout << std::endl;
			//}
			outStream << dataunitio::highQualityCBR(bytes, sliceScalar); // Write output in HQ CBR mode
			outStream << outWrapped;
			if (!outStream) {
				cerr << "Failed to write output file \"" << outFileName << "\"" << endl;
				return EXIT_FAILURE;
			}
		}

		cout << "Encode HQ CBR Done" << endl;
#if 0
		cout << "Please input any kety to exit :";
		cin.get();
#endif

	return EXIT_SUCCESS;
}
