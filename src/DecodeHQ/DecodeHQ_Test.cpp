/***********************************************************************/
/* DecodeHQ.cpp                                                        */
/* Author: Tim Borer                                                   */
/* This version 4th November 2013                                      */
/*                                                                     */
/* Reads compressed transform data in                                  */
/* Decompresses image using VC-2 High Quality profile                  */
/* Writes image data out to a planar file.                             */
/* It is not necessarily complet nor korrect.                          */
/* Copyright (c) BBC 2011-2015 -- For license see the LICENSE file     */
/***********************************************************************/

const char version[] = __DATE__ " @ " __TIME__;
const char summary[] = "Decodes the compressed bytes of a VC-2 High Quality profile to an uncompressed planar file";
const char description[] = "\
This program decodes SMPTE VC-2 HQ profile compressed transform data to regenerate an image sequence.\n\
Its primary output is the decoded image sequence. However it may produce alternative outputs which are:\n\
  1 the wavelet transform of the decoded output (inverse quantised wavelet coefficients)\n\
  2 the quantised wavelet coefficients\n\
  3 the quantisation indices used for each slice\n\
  4 the decoded sequence\n\
Input is just a sequence of compressed bytes.\n\
Output (where appropriate) are in planar format (4:4:4, 4:2:2, 4:2:0 or RGB).\n\
There can be 1 to 4 bytes per sample and the data is left (MSB) justified.\n\
Data is assumed offset binary (which is fine for both YCbCr or RGB).\n\
\n\
Example: DecodeHQ -v -x 1920 -y 1080 -f 4:2:2 -i -l 10 -k LeGall -d 3 -u 1 -a 2 inFileName outFileName";
const char* details[] = { version, summary, description };

#include <cstdlib> //for EXIT_SUCCESS, EXIT_FAILURE, atoi
#include <stdexcept> //For standard logic errors
#include <iostream> //For cin, cout, cerr
#include <string>
#include <fstream>
#include <cstdio> // for perror

#include "DecodeParams.h"
#include "Arrays.h"
#include "Slices.h"
#include "Picture.h"
#include "Frame.h"
#include "Quantisation.h"
#include "WaveletTransform.h"
#include "Utils.h"

using std::cout;
using std::cin;
using std::cerr;
using std::clog;
using std::endl;
using std::string;
using std::filebuf;
using std::streambuf;
using std::ios_base;
using std::istream;
using std::ostream;
using std::fstream;
using std::fill_n;

using arrayio::ioFormat;  // enu
using arrayio::format;    // class
using arrayio::wordWidth;  //class
using arrayio::bitDepth;   //class
using arrayio::offset;    //class

int main(void) {

	int bits = 8;
	const ColourFormat chromaFormat = CF422;

	const string inFileName = "../EncodeHQ_CBR/Haar1D1YUV422Bits10.drc";
	const string outFileName = "Haar1D1YUV422Bits10.ppm";

	float CompressedRate = 2.0;
	const WaveletKernel kernel = LeGall;// Haar1; // {DD97, LeGall, DD137, Haar0, Haar1, Fidelity, Daub97, NullKernel};
	const int waveletDepth = 3;
	const bool verbose = 1;
	const int height = 256;
	const int width = 256;
	const int lumaDepth = bits;
	const int chromaDepth = bits;
	const bool interlaced = 0;
	const bool topFieldFirst = 0;
	const Output output = DECODED;
	int sliceScalar = 1;
	int ySize;
	int xSize;
	int compressedBytes;
	int nbytes;

	if (bits <= 8)
	{
		nbytes = 1;
	}
	else
	{
		nbytes = 2;
	}

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
	default:
		cerr << " Color Format ERROR " << endl;
		return EXIT_FAILURE;
	}



	if (verbose) {
		clog << "input file = " << inFileName << endl;
		clog << "output file = " << outFileName << endl;
	}

	// Open input file or use standard input
	// Input stream is read only binary mode.
	// No point in continuing if can't open input file.
	filebuf inFileBuffer; // For file input. Needs to be defined here to remain in scope
	streambuf *pInBuffer; // Either standard input buffer or a file buffer

	pInBuffer = inFileBuffer.open(inFileName.c_str(), ios_base::in | ios_base::binary);
	if (!pInBuffer) {
		perror((string("Failed to open input file \"") + inFileName + "\"").c_str());
		return EXIT_FAILURE;
	}
	istream inStream(pInBuffer);

	// Open output file or use standard output.
	// Output stream is write only binary mode
	// No point in continuing if can't open output file.
	filebuf outFileBuffer; // For file output. Needs to be defined here to remain in scope
	streambuf *pOutBuffer; // Either standard output buffer or a file buffer
	pOutBuffer = outFileBuffer.open(outFileName.c_str(), ios_base::out | ios_base::binary);
	if (!pOutBuffer) {
		perror((string("Failed to open output file \"") + outFileName + "\"").c_str());
		return EXIT_FAILURE;
	}
	ostream outStream(pOutBuffer);

	if (verbose) {
		clog << "bytes per sample= " << nbytes << endl;
		clog << "luma depth (bits) = " << lumaDepth << endl;
		clog << "chroma depth (bits) = " << chromaDepth << endl;
		clog << "height = " << height << endl;
		clog << "width = " << width << endl;
		clog << "chroma format = " << chromaFormat << endl;
		clog << "interlaced = " << std::boolalpha << interlaced << endl;
		if (interlaced) clog << "top field first = " << std::boolalpha << topFieldFirst << endl;
		clog << "wavelet kernel = " << kernel << endl;
		clog << "wavelet depth = " << waveletDepth << endl;
		clog << "vertical slice size (in units of 2**(wavelet depth)) = " << ySize << endl;
		clog << "horizontal slice size (in units of 2**(wavelet depth)) = " << xSize << endl;
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


		if (verbose) clog << "Determine quantisation indices" << endl;
		const int pictureBytes = (interlaced ? compressedBytes / 2 : compressedBytes);
		// Calculate number of bytes for each slice
		const Array2D bytes = slice_bytes(ySlices, xSlices, pictureBytes, sliceScalar);

		// Construct an container to read the compressed data into.
		const PictureFormat transformFormat(paddedPictureHeight, paddedWidth, chromaFormat);
		Slices inSlices(transformFormat, waveletDepth, ySlices, xSlices);

		// Define picture format (field or frame)
		const PictureFormat picFormat(pictureHeight, width, chromaFormat);

		// Create Frame to hold output data
		const PictureFormat frameFormat(height, width, chromaFormat);
		Frame outFrame(frameFormat, interlaced, topFieldFirst);

		int frame = 0;

		inStream >> sliceio::highQualityCBR(bytes, sliceScalar); // Read input in HQ VBR mode
		inStream >> inSlices; // Read the compressed input picture
								// Check picture was read OK
		if (!inStream) {
			if (frame == 0) {
				cerr << "\rFailed to read the first compressed frame" << endl;
				return EXIT_FAILURE;
			}
			else {
				if (verbose) clog << "\rEnd of input reached after " << frame << " frames     " << endl;
				if (inFileName != "-") inFileBuffer.close();
				if (outFileName != "-") outFileBuffer.close();
				return EXIT_SUCCESS;
			}
		}
		else clog << endl;

		// Reorder quantised coefficients from slice order to transform order
		if (verbose) clog << "Merge slices into full picture" << endl;
		const Picture yuvQCoeffs = merge_blocks(inSlices.yuvSlices);


		// Inverse quantise in transform order
		if (verbose) clog << "Inverse quantise" << endl;
		const Picture yuvTransform = inverse_quantise_transform_np(yuvQCoeffs, inSlices.qIndices, qMatrix);
		
		// Inverse wavelet transform
		if (verbose) clog << "Inverse transform" << endl;
		const Picture outPicture = inverseWaveletTransform(yuvTransform, kernel, waveletDepth, picFormat);

		const Shape2D  restoredSize = { { height, width } };
		Array2D restoredR(restoredSize);
		Array2D restoredG(restoredSize);
		Array2D restoredB(restoredSize);

		Array2D restoredY;
		Array2D restoredU;
		Array2D restoredV;

		const int MaxValue = utils::pow(2, bits) - 1;

		restoredY = clip(outPicture.y(), 0, MaxValue);
		restoredU = clip(outPicture.c1(), 0, MaxValue);
		restoredV = clip(outPicture.c2(), 0, MaxValue);

		int *ULine = (new int[width + 2]) + 1;
		int *VLine = (new int[width + 2]) + 1;
		fill_n(&ULine[-1], width + 2, 0);
		fill_n(&VLine[-1], width + 2, 0);

		const int UVHeight = height + 2;
		const int UVWidth = width + 2;
		const int UVImageSize = UVHeight*UVWidth;
		int *UImage = (new int[UVImageSize]) + UVWidth + 1;
		int *VImage = (new int[UVImageSize]) + UVWidth + 1;
		fill_n(&UImage[-(UVWidth + 1)], UVImageSize, 0);
		fill_n(&VImage[-(UVWidth + 1)], UVImageSize, 0);
		int R, G, B;
		int Y, U, V;

		switch (chromaFormat) {
		case RGB:
			restoredR = clip(outPicture.y(), 0, MaxValue);
			restoredG = clip(outPicture.c1(), 0, MaxValue);
			restoredB = clip(outPicture.c2(), 0, MaxValue);
			break;
		case CF444:
			for (int line = 0; line<height; ++line) {
				for (int pixel = 0; pixel<width; ++pixel) {

					Y = restoredY[line][pixel] - 16;
					U = restoredU[line][pixel] - 128;
					V = restoredV[line][pixel] - 128;

					R = ((298 * Y + 409 * V + 128) >> 8);
					G = ((298 * Y - 100 * U - 208 * V + 128) >> 8);
					B = ((298 * Y + 516 * U + 128) >> 8);

					//Clip RGB Values
					restoredR[line][pixel] = static_cast<int>((R<0) ? 0 : ((R>MaxValue) ? MaxValue : R));
					restoredG[line][pixel] = static_cast<int>((G<0) ? 0 : ((G>MaxValue) ? MaxValue : G));
					restoredB[line][pixel] = static_cast<int>((B<0) ? 0 : ((B>MaxValue) ? MaxValue : B));
				}
			}
			break;
		case CF422:
			for (int line = 0; line<height; ++line) {
				int UVIndex = width*line / 2;
				for (int pixel = 0; pixel<width; pixel += 2) {
					//Copy (sub-sampled) UV components to line buffer.
					ULine[pixel] = restoredU[line][pixel / 2] - 128;
					VLine[pixel] = restoredV[line][pixel / 2] - 128;
				}

				int YIndex = width*line;
				int RGBIndex = 3 * width*line;
				for (int pixel = 0; pixel<width; ++pixel) {

					//Copy Y value and  filter UV values.
					Y = restoredY[line][pixel] - 16;
					U = (ULine[pixel - 1] + 2 * ULine[pixel] + ULine[pixel + 1] + 1) >> 1;
					V = (VLine[pixel - 1] + 2 * VLine[pixel] + VLine[pixel + 1] + 1) >> 1;

					//Matrix YUV to RGB
					R = ((298 * Y + 409 * V + 128) >> 8);
					G = ((298 * Y - 100 * U - 208 * V + 128) >> 8);
					B = ((298 * Y + 516 * U + 128) >> 8);

					//Clip RGB Values
					restoredR[line][pixel] = static_cast<int>((R<0) ? 0 : ((R>MaxValue) ? MaxValue : R));
					restoredG[line][pixel] = static_cast<int>((G<0) ? 0 : ((G>MaxValue) ? MaxValue : G));
					restoredB[line][pixel] = static_cast<int>((B<0) ? 0 : ((B>MaxValue) ? MaxValue : B));
				}
			}
			break;
		case CF420:
			for (int line = 0; line<height; line += 2) {
				for (int pixel = 0; pixel<width; pixel += 2) {
					UImage[line*UVWidth + pixel] = restoredU[line / 2][pixel / 2] - 128;
					VImage[line*UVWidth + pixel] = restoredV[line / 2][pixel / 2] - 128;
				}
			}

			//Vertically interpolate the UV samples
			for (int line = 1; line<height; line += 2) {
				for (int pixel = 0; pixel<width; pixel += 2) {
					UImage[line*UVWidth + pixel] = ((UImage[(line - 1)*UVWidth + pixel] +
						2 * UImage[line*UVWidth + pixel] + UImage[(line + 1)*UVWidth + pixel] + 1) >> 1);
					VImage[line*UVWidth + pixel] = ((VImage[(line - 1)*UVWidth + pixel] +
						2 * VImage[line*UVWidth + pixel] + VImage[(line + 1)*UVWidth + pixel] + 1) >> 1);
				}
			}

			for (int line = 0; line<height; ++line) {
				for (int pixel = 0; pixel<width; ++pixel) {

					//Copy Y value and  filter UV values.
					Y = restoredY[line][pixel] - 16;
					U = (UImage[line*UVWidth + pixel - 1] + 2 * UImage[line*UVWidth + pixel] + UImage[line*UVWidth + pixel + 1] + 1) >> 1;
					V = (VImage[line*UVWidth + pixel - 1] + 2 * VImage[line*UVWidth + pixel] + VImage[line*UVWidth + pixel + 1] + 1) >> 1;

					//Matrix YUV to RGB
					R = ((298 * Y + 409 * V + 128) >> 8);
					G = ((298 * Y - 100 * U - 208 * V + 128) >> 8);
					B = ((298 * Y + 516 * U + 128) >> 8);

					//Clip RGB Values
					restoredR[line][pixel] = static_cast<int>((R<0) ? 0 : ((R>MaxValue) ? MaxValue : R));
					restoredG[line][pixel] = static_cast<int>((G<0) ? 0 : ((G>MaxValue) ? MaxValue : G));
					restoredB[line][pixel] = static_cast<int>((B<0) ? 0 : ((B>MaxValue) ? MaxValue : B));
				}
			}
			break;
		default:
			cerr << "output Image format ERROR" << endl;
			return EXIT_FAILURE;
		}

		delete[](&VLine[-1]);
		delete[](&ULine[-1]);
		delete[](&VImage[-(UVWidth + 1)]);
		delete[](&UImage[-(UVWidth + 1)]);

		if (output == DECODED) {

			string headP6 = "P6";

			outStream << headP6 << endl;
			outStream << height << " " << width << endl;
			outStream << MaxValue << endl;

			outStream << pictureio::wordWidth(nbytes); // Set number of bytes per value in file
			outStream << pictureio::right_justified;
			outStream << pictureio::offset_binary;
			outStream << pictureio::bitDepth(bits, bits); // Set luma and chroma bit depths

														  //Write pixel data line by line
														  //(starting at the botom of the frame because bitmaps are stored upside down!)
			std::streambuf& outbuf = *(outStream.rdbuf());
			int outBufferSize = 3 * width*nbytes;
			unsigned char *outlineBuffer = new unsigned char[outBufferSize];

			unsigned char temp;
			for (int line = 0; line<height; line++) {
				int bufferOffset = 0;
				for (register int pixel = 0; pixel<width; ++pixel) {
					//read RGB values
					R = restoredR[line][pixel];
					G = restoredG[line][pixel];
					B = restoredB[line][pixel];

					//Store RGB values
					if (nbytes == 1)
					{
						outlineBuffer[bufferOffset++] = R;
						outlineBuffer[bufferOffset++] = G;
						outlineBuffer[bufferOffset++] = B;
					}
					else if (nbytes == 2)
					{
						outlineBuffer[bufferOffset++] = R >> 8;
						outlineBuffer[bufferOffset++] = R;
						outlineBuffer[bufferOffset++] = G >> 8;
						outlineBuffer[bufferOffset++] = G;
						outlineBuffer[bufferOffset++] = B >> 8;
						outlineBuffer[bufferOffset++] = B;
					}
					else
					{
						cerr << "Error: nbytes must be 1 or 2" << endl;
						return EXIT_FAILURE;
					}

				} //end pixel loop
				if ((outbuf.sputn(reinterpret_cast<char*>(outlineBuffer), outBufferSize)) < outBufferSize) {
					cerr << "Error: failed to write line " << line << endl;
					return EXIT_FAILURE;
				}
			} //end line loop


			delete[] outlineBuffer;
		}// if (output== DECODE)


		inFileBuffer.close();
		outFileBuffer.close();

		cout << "Decode HD CBR Done " <<endl;
#ifdef _WIN32
		cout << "Please input a character : ";
		cin.get();
#endif

	return EXIT_SUCCESS;
}
