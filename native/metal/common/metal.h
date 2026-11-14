#pragma once
// Metal device, pipelines, buffers and readback for the shared application.

#import <Metal/Metal.h>
#import <QuartzCore/CAMetalLayer.h>

#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

// Must match `View` in shaders.h. Double values are split into float pairs (hi, lo).
struct ViewConstants {
    uint32_t start, end, width, height;
    float cx[2], cy[2], scale[2];
    uint32_t mode, function;
    float density;
    uint32_t sphere;      // 0: complex plane, 1: Riemann sphere
    float radius[2];      // sphere radius in pixels
    float rot[9][2];      // sphere rotation: images of the view x, y, z axes
    uint32_t palette;     // phase colours: 0 hue wheel, 1 GreenPinkTones
    uint32_t count;       // total number of summed items (`end` of them are done when drawing)
};
static_assert(sizeof(ViewConstants) == 144, "Must match the MSL View struct");

// A double as an unevaluated sum of two floats.
inline void splitDouble(double d, float out[2]) {
    out[0] = float(d);
    out[1] = std::isfinite(d) ? float(d - double(out[0])) : 0.0f;
}

class MetalRenderer {
    id<MTLDevice> device = nil;
    id<MTLCommandQueue> queue = nil;
    id<MTLComputePipelineState> computeState[2] = {nil, nil};  // [fp64]
    id<MTLRenderPipelineState> displayState[2] = {nil, nil};
    id<MTLBuffer> acc = nil, mats[2] = {nil, nil};
    std::vector<uint8_t> captured;  // RGB of the last frame drawn with capture
    uint32_t width = 0, height = 0;

    static void check(id<MTLCommandBuffer> commands, const char* operation) {
        [commands waitUntilCompleted];
        if (commands.status == MTLCommandBufferStatusError) {
            char text[512];
            std::snprintf(text, sizeof text, "%s failed: %s", operation,
                          commands.error.localizedDescription.UTF8String);
            fatal(text);
        }
    }
    id<MTLLibrary> compile(const std::string& source, bool fp64) {
        MTLCompileOptions* options = [MTLCompileOptions new];
        // The double-float arithmetic depends on IEEE rounding of every float
        // operation; fast math would reassociate it away.
#if defined(__MAC_15_0) && __MAC_OS_X_VERSION_MAX_ALLOWED >= __MAC_15_0
        if (@available(macOS 15.0, *)) {
            options.mathMode = MTLMathModeSafe;
            options.mathFloatingPointFunctions = MTLMathFloatingPointFunctionsPrecise;
        } else
#endif
        {
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wdeprecated-declarations"
            options.fastMathEnabled = NO;
#pragma clang diagnostic pop
        }
        NSError* error = nil;
        id<MTLLibrary> library = [device newLibraryWithSource:@(source.c_str()) options:options error:&error];
        if (!library) {
            std::fprintf(stderr, "%s\n", error.localizedDescription.UTF8String);
            fatal(fp64 ? "Compiling the fp64 shaders failed" : "Compiling the fp32 shaders failed");
        }
        return library;
    }

public:
    bool ready() const { return acc != nil; }
    id<MTLDevice> metalDevice() const { return device; }

    // sources[0]: fp32 frame, sources[1]: double-float frame (shaderSource in shaders.h)
    std::string initialize(id<MTLDevice> dev, MTLPixelFormat format, const std::string sources[2]) {
        device = dev;
        if (!device) fatal("No Metal device found");
        queue = [device newCommandQueue];
        for (int fp64 = 0; fp64 < 2; fp64++) {
            id<MTLLibrary> library = compile(sources[fp64], fp64);
            NSError* error = nil;
            computeState[fp64] = [device newComputePipelineStateWithFunction:[library newFunctionWithName:@"computeKernel"]
                                                                       error:&error];
            if (!computeState[fp64]) fatal("Create compute pipeline");
            MTLRenderPipelineDescriptor* pd = [MTLRenderPipelineDescriptor new];
            pd.vertexFunction = [library newFunctionWithName:@"vertexMain"];
            pd.fragmentFunction = [library newFunctionWithName:@"fragmentMain"];
            pd.colorAttachments[0].pixelFormat = format;
            displayState[fp64] = [device newRenderPipelineStateWithDescriptor:pd error:&error];
            if (!displayState[fp64]) fatal("Create display pipeline");
        }
        return device.name.UTF8String;
    }
    void resize(uint32_t w, uint32_t h) {
        if (w == width && h == height && acc) return;
        acc = [device newBufferWithLength:size_t(w) * h * 32 options:MTLResourceStorageModeShared];
        if (!acc) fatal("Cannot allocate the accumulation buffer");
        std::memset(acc.contents, 0, acc.length);
        width = w; height = h;
    }
    // f: 8 floats per element; d: 8 doubles per element, uploaded as hi/lo pairs.
    void matrices(const std::vector<float>& f, const std::vector<double>& d) {
        std::vector<float> pairs(d.size() * 2);
        for (size_t i = 0; i < d.size(); i++) splitDouble(d[i], &pairs[2 * i]);
        mats[0] = [device newBufferWithBytes:f.data() length:f.size() * sizeof(float) options:MTLResourceStorageModeShared];
        mats[1] = [device newBufferWithBytes:pairs.data() length:pairs.size() * sizeof(float) options:MTLResourceStorageModeShared];
        if (!mats[0] || !mats[1]) fatal("Cannot allocate the matrix buffers");
    }
    void compute(const ViewConstants& view, bool fp64) {
        id<MTLCommandBuffer> commands = [queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [commands computeCommandEncoder];
        [enc setComputePipelineState:computeState[fp64]];
        [enc setBytes:&view length:sizeof view atIndex:0];
        [enc setBuffer:mats[fp64] offset:0 atIndex:1];
        [enc setBuffer:acc offset:0 atIndex:2];
        [enc dispatchThreads:MTLSizeMake(width, height, 1) threadsPerThreadgroup:MTLSizeMake(16, 16, 1)];
        [enc endEncoding];
        [commands commit];
        check(commands, "GPU summation");
    }
    // Draw into the drawable; with capture, also keep its pixels for pixels().
    void draw(const ViewConstants& view, bool fp64, id<CAMetalDrawable> drawable, bool capture) {
        id<MTLTexture> target = drawable.texture;
        MTLRenderPassDescriptor* pass = [MTLRenderPassDescriptor renderPassDescriptor];
        pass.colorAttachments[0].texture = target;
        pass.colorAttachments[0].loadAction = MTLLoadActionDontCare;
        pass.colorAttachments[0].storeAction = MTLStoreActionStore;
        id<MTLCommandBuffer> commands = [queue commandBuffer];
        id<MTLRenderCommandEncoder> enc = [commands renderCommandEncoderWithDescriptor:pass];
        [enc setRenderPipelineState:displayState[fp64]];
        [enc setFragmentBytes:&view length:sizeof view atIndex:0];
        [enc setFragmentBuffer:acc offset:0 atIndex:1];
        [enc drawPrimitives:MTLPrimitiveTypeTriangle vertexStart:0 vertexCount:3];
        [enc endEncoding];
        id<MTLBuffer> readback = nil;
        const size_t w = target.width, h = target.height, rowBytes = w * 4;
        if (capture) {
            readback = [device newBufferWithLength:rowBytes * h options:MTLResourceStorageModeShared];
            id<MTLBlitCommandEncoder> blit = [commands blitCommandEncoder];
            [blit copyFromTexture:target sourceSlice:0 sourceLevel:0 sourceOrigin:MTLOriginMake(0, 0, 0)
                       sourceSize:MTLSizeMake(w, h, 1) toBuffer:readback destinationOffset:0
              destinationBytesPerRow:rowBytes destinationBytesPerImage:rowBytes * h];
            [blit endEncoding];
        }
        [commands presentDrawable:drawable];
        [commands commit];
        check(commands, "Display");
        if (capture) {
            const uint8_t* bgra = (const uint8_t*)readback.contents;
            captured.resize(w * h * 3);
            for (size_t i = 0; i < w * h; i++) {
                captured[3 * i] = bgra[4 * i + 2];
                captured[3 * i + 1] = bgra[4 * i + 1];
                captured[3 * i + 2] = bgra[4 * i];
            }
        }
    }
    // The two accumulated complex sums at pixel (x, y), as doubles.
    void readPixel(uint32_t x, uint32_t y, double values[4]) {
        const float* a = (const float*)acc.contents + (size_t(y) * width + x) * 8;
        for (int k = 0; k < 4; k++) values[k] = double(a[k]) + double(a[4 + k]);
    }
    const std::vector<uint8_t>& pixels() const { return captured; }
};
