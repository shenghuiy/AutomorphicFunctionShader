#pragma once
// Minimal PNG writer: 8-bit RGB, uncompressed deflate, no dependencies.

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <vector>

static uint32_t crc32(const uint8_t* d, size_t n, uint32_t c = 0xFFFFFFFFu) {
    static uint32_t T[256];
    static bool init = false;
    if (!init) {
        for (uint32_t i = 0; i < 256; i++) {
            uint32_t r = i;
            for (int k = 0; k < 8; k++) r = (r & 1) ? 0xEDB88320u ^ (r >> 1) : r >> 1;
            T[i] = r;
        }
        init = true;
    }
    for (size_t i = 0; i < n; i++) c = T[(c ^ d[i]) & 0xFF] ^ (c >> 8);
    return c;
}

static bool writePNG(const char* path, int w, int h, const std::vector<uint8_t>& rgbTopDown) {
    std::vector<uint8_t> raw;
    raw.reserve(size_t(h) * (1 + 3 * w));
    for (int y = 0; y < h; y++) {
        raw.push_back(0);
        raw.insert(raw.end(), rgbTopDown.begin() + size_t(y) * 3 * w, rgbTopDown.begin() + size_t(y + 1) * 3 * w);
    }
    std::vector<uint8_t> z = {0x78, 0x01};
    uint32_t a = 1, b = 0;
    for (uint8_t v : raw) { a = (a + v) % 65521; b = (b + a) % 65521; }
    for (size_t off = 0; off < raw.size() || off == 0; off += 65535) {
        size_t n = std::min<size_t>(65535, raw.size() - off);
        bool last = off + n >= raw.size();
        z.push_back(last ? 1 : 0);
        z.push_back(uint8_t(n & 0xFF)); z.push_back(uint8_t(n >> 8));
        z.push_back(uint8_t(~n & 0xFF)); z.push_back(uint8_t((~n >> 8) & 0xFF));
        z.insert(z.end(), raw.begin() + off, raw.begin() + off + n);
        if (last) break;
    }
    uint32_t adler = (b << 16) | a;
    for (int s = 24; s >= 0; s -= 8) z.push_back(uint8_t(adler >> s));

    FILE* f = std::fopen(path, "wb");
    if (!f) return false;
    auto be32 = [&](uint32_t v) { uint8_t t[4] = {uint8_t(v >> 24), uint8_t(v >> 16), uint8_t(v >> 8), uint8_t(v)}; std::fwrite(t, 1, 4, f); };
    auto chunk = [&](const char* type, const std::vector<uint8_t>& data) {
        be32(uint32_t(data.size()));
        std::vector<uint8_t> td(type, type + 4);
        td.insert(td.end(), data.begin(), data.end());
        std::fwrite(td.data(), 1, td.size(), f);
        be32(crc32(td.data(), td.size()) ^ 0xFFFFFFFFu);
    };
    const uint8_t sig[8] = {0x89, 'P', 'N', 'G', 0x0D, 0x0A, 0x1A, 0x0A};
    std::fwrite(sig, 1, 8, f);
    std::vector<uint8_t> ihdr = {uint8_t(w >> 24), uint8_t(w >> 16), uint8_t(w >> 8), uint8_t(w),
                                 uint8_t(h >> 24), uint8_t(h >> 16), uint8_t(h >> 8), uint8_t(h), 8, 2, 0, 0, 0};
    chunk("IHDR", ihdr);
    chunk("IDAT", z);
    chunk("IEND", {});
    std::fclose(f);
    return true;
}
