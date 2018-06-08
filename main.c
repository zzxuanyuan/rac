/*
 * Copyright (c) 2016-present, Yann Collet, Facebook, Inc.
 * All rights reserved.
 *
 * This source code is licensed under both the BSD-style license (found in the
 * LICENSE file in the root directory of this source tree) and the GPLv2 (found
 * in the COPYING file in the root directory of this source tree).
 * You may select, at your option, one of the above-listed licenses.
 */


/*-************************************
*  Tuning parameters
**************************************/
#ifndef ZSTDCLI_CLEVEL_DEFAULT
#  define ZSTDCLI_CLEVEL_DEFAULT 3
#endif

#ifndef ZSTDCLI_CLEVEL_MAX
#  define ZSTDCLI_CLEVEL_MAX 19   /* without using --ultra */
#endif



/*-************************************
*  Dependencies
**************************************/
#include "platform.h" /* IS_CONSOLE, PLATFORM_POSIX_VERSION */
#include "util.h"     /* UTIL_HAS_CREATEFILELIST, UTIL_createFileList */
#include <stdio.h>    /* fprintf(), stdin, stdout, stderr */
#include <string.h>   /* strcmp, strlen */
#include <errno.h>    /* errno */
#include "fileio.h"   /* stdinmark, stdoutmark, ZSTD_EXTENSION */
#ifndef ZSTD_NOBENCH
#  include "bench.h"  /* BMK_benchFiles, BMK_SetNbSeconds */
#endif
#ifndef ZSTD_NODICT
#  include "dibio.h"  /* ZDICT_cover_params_t, DiB_trainFromFiles() */
#endif
#define ZSTD_STATIC_LINKING_ONLY   /* ZSTD_maxCLevel */
#include "zstd.h"     /* ZSTD_VERSION_STRING */

#include "lz4.h"
#include "lz4hc.h"
#include <vector>

/*-************************************
*  Constants
**************************************/
#define COMPRESSOR_NAME "zstd command line interface"
#ifndef ZSTD_VERSION
#  define ZSTD_VERSION "v" ZSTD_VERSION_STRING
#endif
#define AUTHOR "Yann Collet"
#define WELCOME_MESSAGE "*** %s %i-bits %s, by %s ***\n", COMPRESSOR_NAME, (int)(sizeof(size_t)*8), ZSTD_VERSION, AUTHOR

#define ZSTD_ZSTDMT "zstdmt"
#define ZSTD_UNZSTD "unzstd"
#define ZSTD_CAT "zstdcat"
#define ZSTD_ZCAT "zcat"
#define ZSTD_GZ "gzip"
#define ZSTD_GUNZIP "gunzip"
#define ZSTD_GZCAT "gzcat"
#define ZSTD_LZMA "lzma"
#define ZSTD_UNLZMA "unlzma"
#define ZSTD_XZ "xz"
#define ZSTD_UNXZ "unxz"
#define ZSTD_LZ4 "lz4"
#define ZSTD_UNLZ4 "unlz4"

#define KB *(1 <<10)
#define MB *(1 <<20)
#define GB *(1U<<30)

#define DISPLAY_LEVEL_DEFAULT 2

static const char*    g_defaultDictName = "dictionary";
static const unsigned g_defaultMaxDictSize = 110 KB;
static const int      g_defaultDictCLevel = 3;
static const unsigned g_defaultSelectivityLevel = 9;
static const unsigned g_defaultMaxWindowLog = 27;
#define OVERLAP_LOG_DEFAULT 9999
#define LDM_PARAM_DEFAULT 9999  /* Default for parameters where 0 is valid */
static U32 g_overlapLog = OVERLAP_LOG_DEFAULT;
static U32 g_ldmHashLog = 0;
static U32 g_ldmMinMatch = 0;
static U32 g_ldmHashEveryLog = LDM_PARAM_DEFAULT;
static U32 g_ldmBucketSizeLog = LDM_PARAM_DEFAULT;


/*-************************************
*  Display Macros
**************************************/
#define DISPLAY(...)         fprintf(g_displayOut, __VA_ARGS__)
#define DISPLAYLEVEL(l, ...) { if (g_displayLevel>=l) { DISPLAY(__VA_ARGS__); } }
static int g_displayLevel = DISPLAY_LEVEL_DEFAULT;   /* 0 : no display,  1: errors,  2 : + result + interaction + warnings,  3 : + progression,  4 : + information */
static FILE* g_displayOut;

struct PerfClock {
    double dictClock;
    double compressClock;
    double decompressClock;

    PerfClock() {
        dictClock = 0.0;
        compressClock = 0.0;
        decompressClock = 0.0;
    }
    ~PerfClock() {}
    void resetClock() {
        dictClock = 0.0;
        compressClock = 0.0;
        decompressClock = 0.0;
    }
    void resetDictClock() {
        dictClock = 0.0;
    }
    void resetCompressClock() {
        compressClock = 0.0;
    }
    void resetDecompressClock() {
        decompressClock = 0.0;
    }
};

struct Parameters {
    const char** filenameTable;
    const char** filedictTable;
    int filenameIdx;
    unsigned maxDictSize;
    int cLevel;
    size_t chunkSize;
    size_t blockSize;
    int blockPerStripe;
    
    Parameters() {
        filenameTable = nullptr;
	filedictTable = nullptr;
	filenameIdx = 0;;
	maxDictSize = 0;
	cLevel = 0;
	chunkSize = 0;
	blockSize = 0;
	blockPerStripe = 0;
    }

    ~Parameters() {}
};

static PerfClock g_perfClock = PerfClock();
static Parameters g_parameters = Parameters();

static ZDICT_cover_params_t defaultCoverParams(void)
{
    ZDICT_cover_params_t params;
    memset(&params, 0, sizeof(params));
    params.d = 8;
    params.steps = 4;
    return params;
}

/* longCommandWArg() :
 *  check if *stringPtr is the same as longCommand.
 *  If yes, @return 1 and advances *stringPtr to the position which immediately follows longCommand.
 * @return 0 and doesn't modify *stringPtr otherwise.
 */
static unsigned longCommandWArg(const char** stringPtr, const char* longCommand)
{
	    size_t const comSize = strlen(longCommand);
	        int const result = !strncmp(*stringPtr, longCommand, comSize);
		    if (result) *stringPtr += comSize;
		        return result;
}

/*! readU32FromChar() :
 *  * @return : unsigned integer value read from input in `char` format.
 *   *  allows and interprets K, KB, KiB, M, MB and MiB suffix.
 *    *  Will also modify `*stringPtr`, advancing it to position where it stopped reading.
 *     *  Note : function will exit() program if digit sequence overflows */
static unsigned readU32FromChar(const char** stringPtr)
{
	const char errorMsg[] = "error: numeric value too large";
	unsigned result = 0;
	while ((**stringPtr >='0') && (**stringPtr <='9')) {
		unsigned const max = (((unsigned)(-1)) / 10) - 1;
		if (result > max) printf(errorMsg);
		result *= 10, result += **stringPtr - '0', (*stringPtr)++ ;
	}
	if ((**stringPtr=='K') || (**stringPtr=='M')) {
		unsigned const maxK = ((unsigned)(-1)) >> 10;
		if (result > maxK) printf(errorMsg);;
		result <<= 10;
		if (**stringPtr=='M') {
			if (result > maxK) printf(errorMsg);
			result <<= 10;
		}
		(*stringPtr)++;  /* skip `K` or `M` */
		if (**stringPtr=='i') (*stringPtr)++;
		if (**stringPtr=='B') (*stringPtr)++;
	}
	return result;
}

static void usage() {
	printf("\n\n");
	printf("usage() :\n");
	printf("run [options] filename\n");
	printf("    -h, help command\n");
	printf("    --block-size, block size\n");
	printf("    --maxdict, dictionary size\n");
	printf("    -e#, compression level\n");
	printf("\n\n");
}

static void concat(char* buf, const char* &x, const size_t xsize, const char* &y, const size_t ysize) {
    buf = (char *)malloc(xsize + ysize);
    memcpy(buf, x, xsize);
    memcpy(buf + xsize, y, ysize);
}

static size_t readFile(const char* filename, char* &buf) {
    FILE* fp = fopen(filename, "rb");
    size_t bufSize = -1;
    if(fp) {
        fseek(fp, 0, SEEK_END);
	bufSize = ftell(fp);
        buf = (char *)malloc(bufSize);
        rewind(fp);
        fread(buf, 1, bufSize, fp);
        fclose(fp);
    } else {
        printf("Error: read file\n");
    }
    return bufSize;
}

static void writeFile(const char* filename, char* &buf, const size_t bufSize) {
    FILE* fp = fopen(filename, "wb");
    if(fp) {
        fwrite(buf, 1, bufSize, fp);
        fclose(fp);
    } else {
        printf("Error: write file\n");
    }
}

struct Header {
    uint32_t    numStripe;
    uint32_t    totalRawSize;
    uint32_t    totalCompressedSize;
};

struct StripeHeader {
    uint32_t    numEntry;
    uint32_t    totalRawSize;
    uint32_t    totalCompressedSize;
};

struct StripeEntry {
    uint32_t    offsetOfCompressedData;
    uint32_t    rawSize;
    uint32_t    compressedSize;
};

struct Entry {
    uint32_t    offsetOfCompressedData;
    uint32_t    rawSize;
    uint32_t    compressedSize;
};

int generateDict(char* &dictBuffer, unsigned maxDictSize, int cLevel, const char** &filedictTable, int filenameIdx, size_t chunkSize) {
    /* generate dictionary */
    ZDICT_cover_params_t coverParams = defaultCoverParams();
    ZDICT_params_t zParams;
    zParams.compressionLevel = cLevel;
    zParams.notificationLevel = 1;
    zParams.dictID = 0;
    coverParams.nbThreads = 1;
    coverParams.zParams = zParams;
    int dictSize = 0;
    double executionTime = 0.0;
    dictSize = DiB_trainFromFilesToBuffer(dictBuffer, maxDictSize, filedictTable, filenameIdx, chunkSize, NULL, &coverParams, 1, executionTime);
    if(executionTime == 0.0) printf("Error: dictionary function return a wrong execution time\n");
    g_perfClock.dictClock += executionTime;
    if(!dictBuffer) printf("dictBuffer is empty, the function was wrong!\n");
    return dictSize;
}

static size_t compressWoutDict(char* &dst, const char* &src, const size_t srcSize, int acceleration, int blockSize, int blockPerStripe) {
    std::vector<Entry> entries;
    size_t dstSize = srcSize * 2;
    char* compressedData = (char*)malloc(dstSize);
    int stripeSize = blockSize * blockPerStripe;
    char* cmpBuf = (char*)malloc(stripeSize * 2);
    size_t offset = 0;
    const char* p = src;
    const char* e = p + srcSize;

    clock_t start_compress = clock();
    while(p < e) {
        size_t dist = (size_t) (e - p);
	size_t len = dist < stripeSize ? dist : stripeSize;
        int cmpBytes = LZ4_compress_fast(p, cmpBuf, len, len*2, acceleration);
	memcpy(compressedData + offset, cmpBuf, cmpBytes);
	entries.emplace_back(Entry { (uint32_t) offset, (uint32_t) len, (uint32_t) cmpBytes });
	offset += cmpBytes;
	p += len;
	if(offset > dstSize) printf("Error : decompression goes beyond output buffer size\n");
    }
    clock_t end_compress = clock();
    double cpu_time_used_compress = ((double) (end_compress - start_compress)) / CLOCKS_PER_SEC;
    g_perfClock.compressClock += cpu_time_used_compress;

    Header header;
    header.numStripe = (uint32_t) entries.size();
    header.totalRawSize = (uint32_t) srcSize;
    header.totalCompressedSize = (uint32_t) offset;
    size_t compressedSize = offset;
    size_t resultSize = sizeof(Header) + sizeof(Entry) * entries.size() + compressedSize;
    char* result = (char*)malloc(resultSize);

    size_t resultOffset = 0;
    memcpy(result, &header, sizeof(Header));
    resultOffset += sizeof(Header);
    memcpy(result + resultOffset, entries.data(), sizeof(Entry) * entries.size());
    resultOffset += sizeof(Entry) * entries.size();
    memcpy(result + resultOffset, compressedData, compressedSize); 
    resultOffset += compressedSize;

    dst = result;
    free(cmpBuf);
    free(compressedData);

    return resultOffset;
}

static size_t compressWoutDictHc(char* &dst, const char* &src, const size_t srcSize, int acceleration, int blockSize, int blockPerStripe) {
    std::vector<Entry> entries;
    size_t dstSize = srcSize * 2;
    char* compressedData = (char*)malloc(dstSize);
    int stripeSize = blockSize * blockPerStripe;
    char* cmpBuf = (char*)malloc(stripeSize * 2);
    size_t offset = 0;
    const char* p = src;
    const char* e = p + srcSize;

    clock_t start_compress = clock();
    while(p < e) {
        size_t dist = (size_t) (e - p);
	size_t len = dist < stripeSize ? dist : stripeSize;
        int cmpBytes = LZ4_compress_HC(p, cmpBuf, len, len*2, acceleration);
	memcpy(compressedData + offset, cmpBuf, cmpBytes);
	entries.emplace_back(Entry { (uint32_t) offset, (uint32_t) len, (uint32_t) cmpBytes });
	offset += cmpBytes;
	p += len;
	if(offset > dstSize) printf("Error : compression goes beyond output buffer size\n");
    }
    clock_t end_compress = clock();
    double cpu_time_used_compress = ((double) (end_compress - start_compress)) / CLOCKS_PER_SEC;
    g_perfClock.compressClock += cpu_time_used_compress;

    Header header;
    header.numStripe = (uint32_t) entries.size();
    header.totalRawSize = (uint32_t) srcSize;
    header.totalCompressedSize = (uint32_t) offset;
    size_t compressedSize = offset;
    size_t resultSize = sizeof(Header) + sizeof(Entry) * entries.size() + compressedSize;
    char* result = (char*)malloc(resultSize);

    size_t resultOffset = 0;
    memcpy(result, &header, sizeof(Header));
    resultOffset += sizeof(Header);
    memcpy(result + resultOffset, entries.data(), sizeof(Entry) * entries.size());
    resultOffset += sizeof(Entry) * entries.size();
    memcpy(result + resultOffset, compressedData, compressedSize); 
    resultOffset += compressedSize;

    dst = result;
    free(cmpBuf);
    free(compressedData);

    return resultOffset;
}

static size_t decompressWoutDict(char* &dst, const char* &src, const size_t srcSize, int blockSize, int blockPerStripe, int stripeIdx, int blockIdx) {

    int stripeSize = blockSize * blockPerStripe;
    const Header* header;
    const StripeEntry* stripes;
    const char* compressedDataTop;

    const char* p = (const char*) src;
    header = (const Header*) p;
    p += sizeof(*header);
    stripes = (const StripeEntry*) p;
    p += sizeof(*stripes) * header->numStripe;
    int indexBound = header->numStripe-1;
    if(stripeIdx > indexBound) printf("Error: stripe index goes beyond boundary\n");

    compressedDataTop = (const char*) p;
    int decompressedSize = 0;
    if(stripeIdx < 0) {

        uint32_t totalRawSize = header->totalRawSize;
        uint32_t totalCompressedSize = header->totalCompressedSize;

        int maxOutputSize = stripeSize * 2;
        char* tempBuf = (char*)malloc(maxOutputSize);
        char* result = (char*)malloc(totalRawSize);
        unsigned offset = 0;

        clock_t start_decompress = clock();
        // Decompress all lines in reverse order for testing.
        for(int iStripe = 0; iStripe < (int) header->numStripe; ++iStripe) {
            const StripeEntry* s = &stripes[iStripe];
            const char* compressedData = compressedDataTop + s->offsetOfCompressedData;
            int tempSize = LZ4_decompress_safe(compressedData, tempBuf, s->compressedSize, maxOutputSize);
            memcpy(result + offset, tempBuf, tempSize);
            offset += tempSize;
        }
        clock_t end_decompress = clock();
        double cpu_time_used_decompress = ((double) (end_decompress - start_decompress)) / CLOCKS_PER_SEC;
        g_perfClock.decompressClock += cpu_time_used_decompress;

        dst = result;
        free(tempBuf);
	decompressedSize = offset;
    } else {
        int maxOutputSize = stripeSize * 2;
        char* tempBuf = (char*)malloc(maxOutputSize);
        
        clock_t start_decompress = clock();
        const StripeEntry* s = &stripes[stripeIdx];
        const char* compressedData = compressedDataTop + s->offsetOfCompressedData;
        int tempSize = LZ4_decompress_safe(compressedData, tempBuf, s->compressedSize, maxOutputSize);
        char* blockBuf = (char*)malloc(blockSize);
        size_t blockPos = blockIdx * blockSize;
        memcpy(blockBuf, tempBuf+blockPos, blockSize);
        clock_t end_decompress = clock();
        double cpu_time_used_decompress = ((double) (end_decompress - start_decompress)) / CLOCKS_PER_SEC;
        g_perfClock.decompressClock += cpu_time_used_decompress;

        dst = blockBuf;
	free(tempBuf);
	decompressedSize = blockSize;
    }
    return decompressedSize;
}

static size_t compressStripeWithDict(LZ4_stream_t* &stream, char* &dst, const char* &src, const size_t srcSize, const size_t blockSize) {
    char* dictBuffer = nullptr;
    int dictSize =  generateDict(dictBuffer, g_parameters.maxDictSize, g_parameters.cLevel, g_parameters.filedictTable, g_parameters.filenameIdx, g_parameters.chunkSize);
    if(dictSize != g_parameters.maxDictSize) printf("Error: dictSize != g_parameters.maxDictSize\n");

    // handle invalid dictionary stripe
    if(dictSize < g_parameters.maxDictSize) {
        StripeHeader stripe;
	stripe.numEntry = 0;
	stripe.totalRawSize = (uint32_t) srcSize;
	stripe.totalCompressedSize = (uint32_t) srcSize;
	size_t resultSize = g_parameters.maxDictSize + sizeof(StripeHeader) + srcSize;
        char* result = (char*)malloc(resultSize);
	size_t resultOffset = g_parameters.maxDictSize;
        memcpy(result + resultOffset, &stripe, sizeof(StripeHeader));
        resultOffset += sizeof(StripeHeader);
        memcpy(result + resultOffset, src, srcSize);
        resultOffset += srcSize;
	dst = result;
	return resultOffset;
    }

    std::vector<Entry> entries;
    size_t dstSize = srcSize * 2;
    char* compressedData = (char*)malloc(dstSize);
    size_t distSize = LZ4_compressBound(blockSize);
    char* cmpBuf = (char*)malloc(distSize);
    size_t offset = 0;
    const char* p = src;
    const char* e = p + srcSize;

    clock_t start_compress = clock();
    while(p < e) {
        size_t dist = (size_t) (e - p);
	size_t len = dist < blockSize ? dist : blockSize;
	LZ4_loadDict(stream, (const char*) dictBuffer, (int) dictSize);
	const size_t cmpBytes = LZ4_compress_fast_continue(stream, p, cmpBuf, (int) len, (int) dstSize, 1);
        memcpy(compressedData + offset, cmpBuf, cmpBytes);
        entries.emplace_back(Entry { (uint32_t) offset, (uint32_t) len, (uint32_t) cmpBytes });
	offset += cmpBytes;
	p += len;
    }
    clock_t end_compress = clock();
    double cpu_time_used_compress = ((double) (end_compress - start_compress)) / CLOCKS_PER_SEC;
    g_perfClock.compressClock += cpu_time_used_compress;

    StripeHeader stripe;
    stripe.numEntry = (uint32_t) entries.size();
    stripe.totalRawSize = (uint32_t) srcSize;
    stripe.totalCompressedSize = (uint32_t) offset;
    size_t compressedStripeSize = (uint32_t) offset;
    size_t resultSize = g_parameters.maxDictSize + sizeof(StripeHeader) + sizeof(Entry) * entries.size() + compressedStripeSize;
    char* result = (char*)malloc(resultSize);
    size_t resultOffset = 0;

    memcpy(result, dictBuffer, dictSize);
    resultOffset += dictSize;
    memcpy(result + resultOffset, &stripe, sizeof(StripeHeader));
    resultOffset += sizeof(StripeHeader);
    memcpy(result + resultOffset, entries.data(), sizeof(Entry) * entries.size());
    resultOffset += sizeof(Entry) * entries.size();
    memcpy(result + resultOffset, compressedData, compressedStripeSize);
    resultOffset += compressedStripeSize;

    dst = result;
    free(dictBuffer);
    free(cmpBuf);
    free(compressedData);
    return resultOffset;
}

static size_t compressWithDict(char* &dst, const char* &src, const size_t srcSize, size_t blockSize, int blockPerStripe) {

    std::vector<StripeEntry> stripes;
    size_t dstSize = srcSize * 2;
    char* compressedData = (char*)malloc(dstSize);
    char* cmpBuf = nullptr;
    size_t offset = 0;
    const char* p = src;
    const char* e = p + srcSize;
    LZ4_stream_t* lz4s;
    lz4s = LZ4_createStream();
    size_t stripeSize = blockSize * blockPerStripe;

    while(p < e) {
        printf("%5.2f%%\r", (p-src) * 100.0 / (e-src));
        size_t dist = (size_t) (e - p);
        size_t len = dist < stripeSize ? dist : stripeSize;

	char* ptr = (char*)malloc(len);
	memcpy(ptr, p, len);
    	writeFile("dictionary.file", ptr, len);

        size_t compressedStripeSize = compressStripeWithDict(lz4s, cmpBuf, p, len, blockSize);
        memcpy(compressedData + offset, cmpBuf, compressedStripeSize);
        stripes.emplace_back(StripeEntry { (uint32_t) offset, (uint32_t) len, (uint32_t) compressedStripeSize });

        offset += compressedStripeSize;
        p += len;
        
	if(offset > dstSize) printf("Error : decompression goes beyond output buffer size\n");

    }

    Header header;
    header.numStripe = (uint32_t) stripes.size();
    header.totalRawSize = (uint32_t) srcSize;
    header.totalCompressedSize = (uint32_t) offset;
    size_t compressedSize = offset;
    size_t resultSize = sizeof(Header) + sizeof(StripeEntry) * stripes.size() + compressedSize;
    char* result = (char*)malloc(resultSize);

    size_t resultOffset = 0;
    memcpy(result, &header, sizeof(Header));
    resultOffset += sizeof(Header);
    memcpy(result + resultOffset, stripes.data(), sizeof(StripeEntry) * stripes.size());
    resultOffset += sizeof(StripeEntry) * stripes.size();
    memcpy(result + resultOffset, compressedData, compressedSize); 
    resultOffset += compressedSize;

    dst = result;
    free(cmpBuf);
    free(compressedData);
    LZ4_freeStream(lz4s);

    return resultOffset;
}

static size_t compressStripeWithDictHc(LZ4_streamHC_t* &stream, char* &dst, const char* &src, const size_t srcSize, size_t blockSize) {


    char* dictBuffer = nullptr;
    int dictSize =  generateDict(dictBuffer, g_parameters.maxDictSize, g_parameters.cLevel, g_parameters.filedictTable, g_parameters.filenameIdx, g_parameters.chunkSize);
    if(dictSize != g_parameters.maxDictSize) printf("Error: dictSize != g_parameters.maxDictSize\n");
    // handle invalid dictionary stripe
    if(dictSize < g_parameters.maxDictSize) {
        StripeHeader stripe;
	stripe.numEntry = 0;
	stripe.totalRawSize = (uint32_t) srcSize;
	stripe.totalCompressedSize = (uint32_t) srcSize;
	size_t resultSize = g_parameters.maxDictSize + sizeof(StripeHeader) + srcSize;
        char* result = (char*)malloc(resultSize);
	size_t resultOffset = g_parameters.maxDictSize;
        memcpy(result + resultOffset, &stripe, sizeof(StripeHeader));
        resultOffset += sizeof(StripeHeader);
        memcpy(result + resultOffset, src, srcSize);
        resultOffset += srcSize;
	dst = result;
	return resultOffset;
    }

    std::vector<Entry> entries;
    size_t dstSize = srcSize * 2;
    char* compressedData = (char*)malloc(dstSize);
    size_t distSize = LZ4_compressBound(blockSize);
    char* cmpBuf = (char*)malloc(distSize);
    size_t offset = 0;
    const char* p = src;
    const char* e = p + srcSize;

    clock_t start_compress = clock();
    while(p < e) {
        size_t dist = (size_t) (e - p);
	size_t len = dist < blockSize ? dist : blockSize;
	LZ4_loadDictHC(stream, (const char*) dictBuffer, (int) dictSize);
	const size_t cmpBytes = LZ4_compress_HC_continue(stream, p, cmpBuf, (int) len, (int) dstSize);
        memcpy(compressedData + offset, cmpBuf, cmpBytes);
        entries.emplace_back(Entry { (uint32_t) offset, (uint32_t) len, (uint32_t) cmpBytes });
	offset += cmpBytes;
	p += len;
    }
    clock_t end_compress = clock();
    double cpu_time_used_compress = ((double) (end_compress - start_compress)) / CLOCKS_PER_SEC;
    g_perfClock.compressClock += cpu_time_used_compress;

    StripeHeader stripe;
    stripe.numEntry = (uint32_t) entries.size();
    stripe.totalRawSize = (uint32_t) srcSize;
    stripe.totalCompressedSize = (uint32_t) offset;
    size_t compressedStripeSize = (uint32_t) offset;
    size_t resultSize = g_parameters.maxDictSize + sizeof(StripeHeader) + sizeof(Entry) * entries.size() + compressedStripeSize;
    char* result = (char*)malloc(resultSize);
    size_t resultOffset = 0;
    memcpy(result + resultOffset, dictBuffer, g_parameters.maxDictSize);
    resultOffset += g_parameters.maxDictSize;
    memcpy(result + resultOffset, &stripe, sizeof(StripeHeader));
    resultOffset += sizeof(StripeHeader);
    memcpy(result + resultOffset, entries.data(), sizeof(Entry) * entries.size());
    resultOffset += sizeof(Entry) * entries.size();
    memcpy(result + resultOffset, compressedData, compressedStripeSize);
    resultOffset += compressedStripeSize;

    dst = result;
    free(dictBuffer);
    free(cmpBuf);
    free(compressedData);
    return resultOffset;
}

static size_t compressWithDictHc(char* &dst, const char* &src, const size_t srcSize, size_t blockSize, int blockPerStripe) {

    std::vector<StripeEntry> stripes;
    size_t dstSize = srcSize * 2;
    char* compressedData = (char*)malloc(dstSize);
    char* cmpBuf = nullptr;
    size_t offset = 0;
    const char* p = src;
    const char* e = p + srcSize;
    LZ4_streamHC_t* lz4s;
    lz4s = LZ4_createStreamHC();
    size_t stripeSize = blockSize * blockPerStripe;

    while(p < e) {
        printf("%5.2f%%\r", (p-src) * 100.0 / (e-src));
        size_t dist = (size_t) (e - p);
        size_t len = dist < stripeSize ? dist : stripeSize;

	char* ptr = (char*)malloc(len);
	memcpy(ptr, p, len);
    	writeFile("dictionary.file", ptr, len);

	size_t compressedStripeSize = compressStripeWithDictHc(lz4s, cmpBuf, p, len, blockSize);
        memcpy(compressedData + offset, cmpBuf, compressedStripeSize);
        stripes.emplace_back(StripeEntry { (uint32_t) offset, (uint32_t) len, (uint32_t) compressedStripeSize });

        offset += compressedStripeSize;
        p += len;
        
	if(offset > dstSize) printf("Error : decompression goes beyond output buffer size\n");

    }

    Header header;
    header.numStripe = (uint32_t) stripes.size();
    header.totalRawSize = (uint32_t) srcSize;
    header.totalCompressedSize = (uint32_t) offset;
    size_t compressedSize = offset;
    size_t resultSize = sizeof(Header) + sizeof(StripeEntry) * stripes.size() + compressedSize;
    char* result = (char*)malloc(resultSize);

    size_t resultOffset = 0;
    memcpy(result, &header, sizeof(Header));
    resultOffset += sizeof(Header);
    memcpy(result + resultOffset, stripes.data(), sizeof(StripeEntry) * stripes.size());
    resultOffset += sizeof(StripeEntry) * stripes.size();
    memcpy(result + resultOffset, compressedData, compressedSize); 
    resultOffset += compressedSize;

    dst = result;
    free(cmpBuf);
    free(compressedData);
    LZ4_freeStreamHC(lz4s);

    return resultOffset;
}

static size_t decompressBlockWithDict(char* &dst, const char* &src, const size_t srcSize, const size_t blockSize, const int blockIdx) {
    const StripeHeader* stripe;
    const Entry* entries;
    const char* compressedDataTop;

    const char* p = (const char*) src;
    char* dict = (char*)malloc(g_parameters.maxDictSize);
    memcpy(dict, p, g_parameters.maxDictSize);
    p += g_parameters.maxDictSize;
    stripe = (const StripeHeader*) p;
    p += sizeof(*stripe);

    // handle invalid dict stripe
    if(stripe->numEntry == 0) {
        if(stripe->totalRawSize != stripe->totalCompressedSize) printf("stripe raw size and compressed size are mismatched\n");
	clock_t start_decompress = clock();
        char* result = (char*)malloc(blockSize);
	size_t blockOffset = blockIdx * blockSize;
        memcpy(result, p + blockOffset, blockSize);
        dst = result;
	clock_t end_decompress = clock();
	double cpu_time_used_decompress = ((double) (end_decompress - start_decompress)) / CLOCKS_PER_SEC;
	g_perfClock.decompressClock += cpu_time_used_decompress;
        return blockSize;
    }

    entries = (const Entry*) p;
    p += sizeof(*stripe) * stripe->numEntry;
    if(blockIdx > stripe->numEntry-1) printf("Error: block index goes beyond boundary\n");
    compressedDataTop = (const char*) p;

    int maxOutputSize = blockSize * 2;
    char* tempBuf = (char*)malloc(maxOutputSize);

    clock_t start_decompress = clock();
    // Decompress single entry.
    const Entry* e = &entries[blockIdx];
    const char* compressedData = compressedDataTop + e->offsetOfCompressedData;
    const int decompressedSize = LZ4_decompress_safe_usingDict(
          (const char*) compressedData
        , tempBuf
        , e->compressedSize
        , maxOutputSize
        , dict
	, (int) g_parameters.maxDictSize);

    if(decompressedSize != blockSize) printf("Error: decompressed size != block size\n");
    char* result = (char*)malloc(blockSize);
    memcpy(result, tempBuf, decompressedSize);
    clock_t end_decompress = clock();
    double cpu_time_used_decompress = ((double) (end_decompress - start_decompress)) / CLOCKS_PER_SEC;
    g_perfClock.decompressClock += cpu_time_used_decompress;

    dst = result;
    free(dict);
    free(tempBuf);
    return decompressedSize;
}

static size_t decompressStripeWithDict(char* &dst, const char* &src, const size_t srcSize, const size_t blockSize) {
    const StripeHeader* stripe;
    const Entry* entries;
    const char* compressedDataTop;

    const char* p = (const char*) src;
    char* dict = (char*)malloc(g_parameters.maxDictSize);
    memcpy(dict, p, g_parameters.maxDictSize);
    p += g_parameters.maxDictSize;
    stripe = (const StripeHeader*) p;
    p += sizeof(*stripe);

    // handle invalid dict stripe
    if(stripe->numEntry == 0) {
        if(stripe->totalRawSize != stripe->totalCompressedSize) printf("stripe raw size and compressed size are mismatched\n");
        char* result = (char*)malloc(stripe->totalRawSize);
	memcpy(result, p, stripe->totalRawSize);
	dst = result;
	return stripe->totalRawSize;
    }

    entries = (const Entry*) p;
    p += sizeof(*stripe) * stripe->numEntry;

    compressedDataTop = (const char*) p;

    uint32_t rawStripeSize = stripe->totalRawSize;
    uint32_t compressedStripeSize = stripe->totalCompressedSize;

    int maxOutputSize = blockSize * 2;
    char* tempBuf = (char*)malloc(maxOutputSize);
    char* result = (char*)malloc(rawStripeSize);
    unsigned offset = 0;

    clock_t start_decompress = clock();
    // Decompress all lines in reverse order for testing.
    for(int iEntry = 0; iEntry < (int) stripe->numEntry; ++iEntry) {
        const Entry* e = &entries[iEntry];
        const char* compressedData = compressedDataTop + e->offsetOfCompressedData;
        const int decompressedSize = LZ4_decompress_safe_usingDict(
              (const char*) compressedData
            , tempBuf
            , e->compressedSize
            , maxOutputSize
            , dict
            , (int) g_parameters.maxDictSize);

	memcpy(result + offset, tempBuf, decompressedSize);
	offset += decompressedSize;
    }
    clock_t end_decompress = clock();
    double cpu_time_used_decompress = ((double) (end_decompress - start_decompress)) / CLOCKS_PER_SEC;
    g_perfClock.decompressClock += cpu_time_used_decompress;

    dst = result;
    free(dict);
    free(tempBuf);
    return offset;
}

static size_t decompressWithDict(char* &dst, const char* &src, const size_t srcSize, const size_t blockSize, const size_t blockPerStripe, int stripeIdx, int blockIdx) {
    const size_t stripeSize = blockSize * blockPerStripe;

    const Header* header;
    const StripeEntry* stripes;
    const char* compressedDataTop;

    const char* p = (const char*) src;
    header = (const Header*) p;
    p += sizeof(*header);
    stripes = (const StripeEntry*) p;
    p += sizeof(*stripes) * header->numStripe;
    int indexBound = header->numStripe-1;
    if(stripeIdx > indexBound) printf("Error: stripe index goes beyond boundary\n");

    compressedDataTop = (const char*) p;
    int decompressedSize = 0;
    if(stripeIdx < 0) {
        uint32_t totalRawSize = header->totalRawSize;
        uint32_t totalCompressedSize = header->totalCompressedSize;
        int maxOutputSize = stripeSize*2;
        char* tempBuf = nullptr;
        char* result = (char*)malloc(totalRawSize);

        unsigned offset = 0;
        // Decompress all lines in reverse order for testing.
        for(int iStripe = 0; iStripe < (int) header->numStripe; ++iStripe) {
            const StripeEntry* s = &stripes[iStripe];
            const char* compressedData = compressedDataTop + s->offsetOfCompressedData;
            int tempSize = decompressStripeWithDict(tempBuf, compressedData, s->compressedSize, blockSize);
            memcpy(result + offset, tempBuf, tempSize);
            offset += tempSize;
        }
	free(tempBuf);
	dst = result;
	decompressedSize = offset;
    } else {
        char* tempBuf = nullptr;
        const StripeEntry* s = &stripes[stripeIdx];
        const char* compressedData = compressedDataTop + s->offsetOfCompressedData;
        decompressedSize = decompressBlockWithDict(tempBuf, compressedData, s->compressedSize, blockSize, blockIdx);
        dst = tempBuf;
    }
    return decompressedSize;
}

static int generateStripeIdx(int numStripe) {
    return rand() % numStripe;
}

static int generateBlockIdx(int numBlock) {
    return rand() % numBlock;
}
 
void testLz4(const char** &filenameTable, unsigned filenameIdx, int cLevel, int blockSize, int blockPerStripe)
{
    /* compress data */
    char* compressBuf = nullptr;
    char* srcBuf = nullptr;
    if(filenameIdx > 1) printf("too many file names\n");
    const char* filename = filenameTable[0];
    int srcSize = readFile(filename, srcBuf);
    const char* srcBufConst = srcBuf;
    int compressSize = -1;
    int dstSize = srcSize * 2;
    if(cLevel >= ZSTDCLI_CLEVEL_DEFAULT) {
        compressSize = compressWoutDictHc(compressBuf, srcBufConst, srcSize, cLevel, blockSize, blockPerStripe);
    } else {
        compressSize = compressWoutDict(compressBuf, srcBufConst, srcSize, 1, blockSize, blockPerStripe);
    }
    writeFile("compressed.file.wout.dict", compressBuf, compressSize);
    printf("compression without dictionary took %f seconds to execute \n", g_perfClock.compressClock);

    /* decompress data */
    char* decompressBuf = nullptr;
    const char* compressBufConst = compressBuf;
    int dstCapacity = srcSize * 2;
    int decompressSize = decompressWoutDict(decompressBuf, compressBufConst, compressSize, blockSize, blockPerStripe, -1, -1);
    writeFile("decompress.file.wout.dict", decompressBuf, decompressSize);
    free(decompressBuf);
    printf("sequential decompression without dictionary took %f seconds to execute \n", g_perfClock.decompressClock);

    g_perfClock.resetDecompressClock();
    int totalDecompressedSize = 0;
    /* randomly decompress data */
    size_t stripeSize = blockSize * blockPerStripe;
    int numStripe = srcSize/stripeSize;
    int numBlock = blockPerStripe;
    int numIteration = srcSize/blockSize;
    double cpu_time_used_random = 0.0;
    for(int iteration = 0; iteration < numIteration; ++iteration) {
        int stripeIdx = generateStripeIdx(numStripe);
        int blockIdx = generateBlockIdx(numBlock);
        decompressBuf = nullptr;
        compressBufConst = compressBuf;
        decompressSize = decompressWoutDict(decompressBuf, compressBufConst, compressSize, blockSize, blockPerStripe, stripeIdx, blockIdx);
	totalDecompressedSize += decompressSize;

	/* verify block data */
        size_t pos = stripeIdx * stripeSize + blockIdx * blockSize;
	if(memcmp(decompressBuf, srcBuf+pos, blockSize) != 0) printf("Error: block data is invalid!\n");
	free(decompressBuf);
    }
    printf("random decompression without dictionary random took %f seconds to execute, total decompressed data %d \n", g_perfClock.decompressClock, totalDecompressedSize);
    free(srcBuf);
    free(compressBuf);
}

   
void testRac(const char** &filenameTable, unsigned filenameIdx, int cLevel, unsigned maxDictSize, size_t blockSize, int blockPerStripe) {

    /* compress data */
    char* compressBuf = nullptr;
    char* srcBuf = nullptr;
    if(filenameIdx > 1) printf("too many file names\n");
    const char* filename = filenameTable[0];
    int srcSize = readFile(filename, srcBuf);
    const char* srcBufConst = srcBuf;
    int compressSize = -1;
    if(cLevel >= ZSTDCLI_CLEVEL_DEFAULT) {
        compressSize = compressWithDictHc(compressBuf, srcBufConst, srcSize, blockSize, blockPerStripe);
    } else {
        compressSize = compressWithDict(compressBuf, srcBufConst, srcSize, blockSize, blockPerStripe);
    }
    writeFile("compressed.file.with.dict", compressBuf, compressSize);
    printf("dictionary took %f seconds to execute \n", g_perfClock.dictClock);
    printf("compression with dictionary took %f seconds to execute \n", g_perfClock.compressClock);

    /* sequentially decompress data */
    char* decompressBuf = nullptr;
    const char* compressBufConst = compressBuf;
    int decompressSize = decompressWithDict(decompressBuf, compressBufConst, compressSize, blockSize, blockPerStripe, -1, -1);
    writeFile("decompress.file.with.dict", decompressBuf, decompressSize);
    free(decompressBuf);
    printf("sequential decompression with dictionary took %f seconds to execute \n", g_perfClock.decompressClock);

    g_perfClock.resetDecompressClock();
    int totalDecompressedSize = 0;
    /* randomly decompress data */
    size_t stripeSize = blockSize * blockPerStripe;
    int numStripe = srcSize/stripeSize;
    int numBlock = blockPerStripe;
    int numIteration = srcSize/blockSize;
    for(int iteration = 0; iteration < numIteration; ++iteration) {
        int stripeIdx = generateStripeIdx(numStripe);
        int blockIdx = generateBlockIdx(numBlock);
        decompressBuf = nullptr;
        compressBufConst = compressBuf;
        decompressSize = decompressWithDict(decompressBuf, compressBufConst, compressSize, blockSize, blockPerStripe, stripeIdx, blockIdx);
	totalDecompressedSize += decompressSize;
        
	/* verify block data */
        size_t pos = stripeIdx * stripeSize + blockIdx * blockSize;
	if(memcmp(decompressBuf, srcBuf+pos, blockSize) != 0) printf("Error: block data is invalid!\n");
	free(decompressBuf);
    }
    printf("random decompression with dictionary took %f seconds to execute, total decompressed size %d \n", g_perfClock.decompressClock, totalDecompressedSize);
    free(srcBuf);
    free(compressBuf);
}

int main(int argCount, const char* argv[])
{
    const char** filenameTable = (const char**)malloc(sizeof(const char*));   /* argCount >= 1 */
    unsigned filenameIdx = 0;
//    const char* str = (const char*)malloc(4);
//    str = "big";
//    filenameTable[0] = str;
    int argNb;
    int lastCommand = 0;
    int nextArgumentsAreFiles = 0;
    int nextArgumentIsMaxDict = 0;
    int cLevel = 3;

    size_t blockSize = 0;
    size_t chunkSize = 0;
    unsigned maxDictSize = 102400;
    unsigned blockPerStripe = 8;
    
    /* command switches */
    for (argNb=1; argNb<argCount; argNb++) {
        const char* argument = argv[argNb];
        if(!argument) continue;   /* Protection if argument empty */

	if (nextArgumentsAreFiles==0) {

            /* Decode commands (note : aggregated commands are allowed) */
            if (argument[0]=='-') {

                if (argument[1]=='-') {
                    /* long commands (--long-word) */
                    if (!strcmp(argument, "--")) { nextArgumentsAreFiles=1; continue; }   /* only file names allowed from now on */
                    if (!strcmp(argument, "--maxdict")) { nextArgumentIsMaxDict=1; lastCommand=1; continue; }  /* kept available for compatibility with old syntax ; will be removed one day */

                    /* long commands with arguments */
                    if (longCommandWArg(&argument, "--dict-chunk-size=")) { chunkSize = readU32FromChar(&argument); continue; }
                    if (longCommandWArg(&argument, "--block-size=")) { blockSize = readU32FromChar(&argument); continue; }
                    if (longCommandWArg(&argument, "--maxdict=")) { maxDictSize = readU32FromChar(&argument); continue; }
                    if (longCommandWArg(&argument, "--blocks-per-stripe=")) { blockPerStripe = readU32FromChar(&argument); continue; }
               }

                argument++;
                while (argument[0]!=0) {
                    if (lastCommand) {
                        DISPLAY("error : command must be followed by argument \n");
                        return 1;
                    }
                    /* compression Level */
                    if ((*argument>='0') && (*argument<='9')) {
                        cLevel = readU32FromChar(&argument);
                        continue;
                    }

                    switch(argument[0])
                    {
                        /* Display help */
                    case 'H':
                    case 'h': g_displayOut=stdout; usage(); return 0;

                        /* range bench (benchmark only) */
                    case 'e':
                        /* compression Level */
                        argument++;
                        cLevel = readU32FromChar(&argument);
                        break;

                        /* cut input into blocks (benchmark only) */
                    case 'B':
                        argument++;
                        blockSize = readU32FromChar(&argument);
                        break;

                        /* unknown command */
                    default : usage(); return 1;
                    }
                }
                continue;
            }   /* if (argument[0]=='-') */

            if (nextArgumentIsMaxDict) {  /* kept available for compatibility with old syntax ; will be removed one day */
                nextArgumentIsMaxDict = 0;
                lastCommand = 0;
                maxDictSize = readU32FromChar(&argument);
                continue;
            }

        }   /* if (nextArgumentIsAFile==0) */

        /* add filename to list */
        filenameTable[filenameIdx++] = argument;
    }

    const char** filedictTable = (const char**)malloc(sizeof(const char*));   /* argCount >= 1 */
    g_parameters.filenameTable = filenameTable;
    const char* filename = "dictionary.file";
    filedictTable[0] = filename;
    g_parameters.filedictTable = filedictTable;
    g_parameters.filenameIdx = filenameIdx;
    g_parameters.cLevel = cLevel;
    g_parameters.maxDictSize = maxDictSize;
    g_parameters.blockSize = blockSize;
    g_parameters.chunkSize = chunkSize;
    g_parameters.blockPerStripe = blockPerStripe;

    g_perfClock.resetClock();
    testRac(filenameTable, filenameIdx, cLevel, maxDictSize, blockSize, blockPerStripe);
    g_perfClock.resetClock();
    testLz4(filenameTable, filenameIdx, cLevel, blockSize, blockPerStripe);
    return 0;
}
