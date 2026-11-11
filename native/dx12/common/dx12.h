#pragma once
#include <d3d12.h>
#include <dxgi1_6.h>
#include <d3dcompiler.h>
#include <wrl/client.h>
#include <cstring>
#include <string>

using Microsoft::WRL::ComPtr;

static void dxCheck(HRESULT hr, const char* operation) {
    if (FAILED(hr)) {
        char text[512];
        std::snprintf(text, sizeof text, "%s failed (HRESULT 0x%08lx)", operation, (unsigned long)hr);
        fatal(text);
    }
}

struct ViewConstants {
    UINT start, end, width, height;
    double cx, cy, scale;
    UINT mode, function;
    float density;
    UINT sphere;          // 0: complex plane, 1: Riemann sphere
    double radius;        // sphere radius in pixels
    double rot[3][3];     // sphere rotation: images of the view x, y, z axes
    UINT palette = 0;     // phase colours: 0 hue wheel, 1 GreenPinkTones
    UINT count = 0;       // total number of summed items (uEnd of them are done when drawing)
};
static_assert(sizeof(ViewConstants) == 144, "Must match the HLSL constant buffer");

class DxRenderer {
    ComPtr<ID3D12Device> device;
    ComPtr<ID3D12CommandQueue> queue;
    ComPtr<ID3D12CommandAllocator> allocator;
    ComPtr<ID3D12GraphicsCommandList> commands;
    ComPtr<ID3D12Fence> fence;
    HANDLE fenceEvent = nullptr;
    UINT64 fenceValue = 0;
    ComPtr<IDXGISwapChain3> swap;
    ComPtr<ID3D12DescriptorHeap> rtvHeap;
    ComPtr<ID3D12Resource> back[2], acc, matF, matD, cursorReadback;
    ComPtr<ID3D12RootSignature> root;
    ComPtr<ID3D12PipelineState> computeF, computeD, display;
    UINT rtvStride = 0, width = 0, height = 0;
    UINT lastBack = 0;
    std::string source;  // full HLSL: common frame + the group's term/value functions

    ComPtr<ID3D12Resource> buffer(UINT64 bytes, D3D12_HEAP_TYPE type,
                                 D3D12_RESOURCE_STATES state, D3D12_RESOURCE_FLAGS flags = D3D12_RESOURCE_FLAG_NONE) {
        D3D12_HEAP_PROPERTIES heap = {};
        heap.Type = type;
        D3D12_RESOURCE_DESC desc = {};
        desc.Dimension = D3D12_RESOURCE_DIMENSION_BUFFER;
        desc.Width = bytes; desc.Height = 1; desc.DepthOrArraySize = 1;
        desc.MipLevels = 1; desc.SampleDesc.Count = 1;
        desc.Layout = D3D12_TEXTURE_LAYOUT_ROW_MAJOR; desc.Flags = flags;
        ComPtr<ID3D12Resource> result;
        dxCheck(device->CreateCommittedResource(&heap, D3D12_HEAP_FLAG_NONE, &desc, state,
                nullptr, IID_PPV_ARGS(&result)), "Create buffer");
        return result;
    }
    void begin() {
        dxCheck(allocator->Reset(), "Reset command allocator");
        dxCheck(commands->Reset(allocator.Get(), nullptr), "Reset command list");
    }
    void finish() {
        dxCheck(commands->Close(), "Close command list");
        ID3D12CommandList* lists[] = {commands.Get()};
        queue->ExecuteCommandLists(1, lists);
        wait();
    }
    void wait() {
        const UINT64 value = ++fenceValue;
        dxCheck(queue->Signal(fence.Get(), value), "Signal GPU fence");
        if (fence->GetCompletedValue() < value) {
            dxCheck(fence->SetEventOnCompletion(value, fenceEvent), "Wait for GPU fence");
            if (WaitForSingleObject(fenceEvent, 30000) != WAIT_OBJECT_0)
                fatal("GPU did not complete work within 30 seconds");
        }
        dxCheck(device->GetDeviceRemovedReason(), "GPU execution");
    }
    void transition(ID3D12Resource* resource, D3D12_RESOURCE_STATES before, D3D12_RESOURCE_STATES after) {
        D3D12_RESOURCE_BARRIER barrier = {};
        barrier.Type = D3D12_RESOURCE_BARRIER_TYPE_TRANSITION;
        barrier.Transition = {resource, D3D12_RESOURCE_BARRIER_ALL_SUBRESOURCES, before, after};
        commands->ResourceBarrier(1, &barrier);
    }
    ComPtr<ID3DBlob> shader(const char* entry, const char* target, bool fp64 = false) {
        D3D_SHADER_MACRO macros[] = {{"COMPUTE","1"}, {"VEC2",fp64 ? "double2" : "float2"},
            {"VEC4",fp64 ? "double4" : "float4"}, {fp64 ? "FP64" : nullptr,"1"}, {nullptr,nullptr}};
        ComPtr<ID3DBlob> code, errors;
        HRESULT hr = D3DCompile(source.data(), source.size(), "plot.hlsl",
            entry[0] == 'c' ? macros : nullptr, nullptr, entry, target,
            D3DCOMPILE_ENABLE_STRICTNESS | D3DCOMPILE_OPTIMIZATION_LEVEL3 | D3DCOMPILE_IEEE_STRICTNESS, 0, &code, &errors);
        if (errors) std::fprintf(stderr, "%s\n", (char*)errors->GetBufferPointer());
        dxCheck(hr, entry);
        return code;
    }
    static D3D12_SHADER_BYTECODE bytecode(ID3DBlob* code) {
        return {code->GetBufferPointer(), code->GetBufferSize()};
    }
    ComPtr<ID3D12Resource> upload(const void* data, size_t bytes) {
        auto staging = buffer(bytes, D3D12_HEAP_TYPE_UPLOAD, D3D12_RESOURCE_STATE_GENERIC_READ);
        void* mapped = nullptr;
        D3D12_RANGE empty = {0,0};
        dxCheck(staging->Map(0, &empty, &mapped), "Map upload buffer");
        std::memcpy(mapped, data, bytes);
        staging->Unmap(0, nullptr);
        auto result = buffer(bytes, D3D12_HEAP_TYPE_DEFAULT, D3D12_RESOURCE_STATE_COPY_DEST);
        begin();
        commands->CopyBufferRegion(result.Get(), 0, staging.Get(), 0, bytes);
        transition(result.Get(), D3D12_RESOURCE_STATE_COPY_DEST, D3D12_RESOURCE_STATE_NON_PIXEL_SHADER_RESOURCE);
        finish();
        return result;
    }

public:
    ~DxRenderer() { if (fenceEvent) CloseHandle(fenceEvent); }
    bool ready() const { return acc != nullptr; }
    std::string initialize(HWND hwnd, UINT w, UINT h, std::string hlsl) {
        source = std::move(hlsl);        // Enable validation when the optional Windows Graphics Tools component is installed.
        ComPtr<ID3D12Debug> debug;
        if (SUCCEEDED(D3D12GetDebugInterface(IID_PPV_ARGS(&debug)))) debug->EnableDebugLayer();
        ComPtr<IDXGIFactory6> factory;
        dxCheck(CreateDXGIFactory2(0, IID_PPV_ARGS(&factory)), "Create DXGI factory");
        DXGI_ADAPTER_DESC1 selected = {};
        for (UINT i=0; ; ++i) {
            ComPtr<IDXGIAdapter1> adapter;
            HRESULT hr=factory->EnumAdapterByGpuPreference(i, DXGI_GPU_PREFERENCE_HIGH_PERFORMANCE, IID_PPV_ARGS(&adapter));
            if (hr == DXGI_ERROR_NOT_FOUND) break;
            dxCheck(hr, "Enumerate adapters");
            DXGI_ADAPTER_DESC1 desc;
            dxCheck(adapter->GetDesc1(&desc), "Get adapter description");
            if (desc.Flags & DXGI_ADAPTER_FLAG_SOFTWARE) continue;
            ComPtr<ID3D12Device> candidate;
            if (FAILED(D3D12CreateDevice(adapter.Get(), D3D_FEATURE_LEVEL_11_0, IID_PPV_ARGS(&candidate)))) continue;
            D3D12_FEATURE_DATA_D3D12_OPTIONS options = {};
            if (FAILED(candidate->CheckFeatureSupport(D3D12_FEATURE_D3D12_OPTIONS, &options, sizeof options)) ||
                !options.DoublePrecisionFloatShaderOps) continue;
            device=candidate; selected=desc; break;
        }
        if (!device) fatal("A DirectX 12 GPU with double-precision shader support is required");
        char gpu[512];
        WideCharToMultiByte(CP_UTF8,0,selected.Description,-1,gpu,sizeof gpu,nullptr,nullptr);
        D3D12_COMMAND_QUEUE_DESC q = {};
        q.Type = D3D12_COMMAND_LIST_TYPE_DIRECT;
        dxCheck(device->CreateCommandQueue(&q, IID_PPV_ARGS(&queue)), "Create command queue");
        dxCheck(device->CreateCommandAllocator(q.Type, IID_PPV_ARGS(&allocator)), "Create allocator");
        dxCheck(device->CreateCommandList(0,q.Type,allocator.Get(),nullptr,IID_PPV_ARGS(&commands)), "Create command list");
        dxCheck(commands->Close(), "Close initial list");
        dxCheck(device->CreateFence(0,D3D12_FENCE_FLAG_NONE,IID_PPV_ARGS(&fence)), "Create fence");
        fenceEvent=CreateEventW(nullptr,FALSE,FALSE,nullptr);
        if (!fenceEvent) fatal("Cannot create GPU fence event");

        DXGI_SWAP_CHAIN_DESC1 sc = {};
        sc.Width=w; sc.Height=h; sc.Format=DXGI_FORMAT_R8G8B8A8_UNORM;
        sc.SampleDesc.Count=1; sc.BufferUsage=DXGI_USAGE_RENDER_TARGET_OUTPUT;
        sc.BufferCount=2; sc.SwapEffect=DXGI_SWAP_EFFECT_FLIP_DISCARD;
        ComPtr<IDXGISwapChain1> chain;
        dxCheck(factory->CreateSwapChainForHwnd(queue.Get(),hwnd,&sc,nullptr,nullptr,&chain), "Create swap chain");
        dxCheck(chain.As(&swap), "Get swap chain");
        dxCheck(factory->MakeWindowAssociation(hwnd,DXGI_MWA_NO_ALT_ENTER), "Configure window association");
        D3D12_DESCRIPTOR_HEAP_DESC hd = {};
        hd.Type=D3D12_DESCRIPTOR_HEAP_TYPE_RTV; hd.NumDescriptors=2;
        dxCheck(device->CreateDescriptorHeap(&hd,IID_PPV_ARGS(&rtvHeap)), "Create render target descriptors");
        rtvStride=device->GetDescriptorHandleIncrementSize(hd.Type);

        D3D12_ROOT_PARAMETER params[4] = {};
        params[0].ParameterType=D3D12_ROOT_PARAMETER_TYPE_32BIT_CONSTANTS;
        params[0].Constants={0,0,sizeof(ViewConstants)/4};
        params[1].ParameterType=D3D12_ROOT_PARAMETER_TYPE_SRV; params[1].Descriptor={0,0};
        params[2].ParameterType=D3D12_ROOT_PARAMETER_TYPE_UAV; params[2].Descriptor={0,0};
        params[3].ParameterType=D3D12_ROOT_PARAMETER_TYPE_SRV; params[3].Descriptor={1,0};
        D3D12_ROOT_SIGNATURE_DESC rd = {};
        rd.NumParameters=4; rd.pParameters=params;
        ComPtr<ID3DBlob> serialized, errors;
        dxCheck(D3D12SerializeRootSignature(&rd,D3D_ROOT_SIGNATURE_VERSION_1,&serialized,&errors), "Serialize root signature");
        dxCheck(device->CreateRootSignature(0,serialized->GetBufferPointer(),serialized->GetBufferSize(),IID_PPV_ARGS(&root)), "Create root signature");
        for (bool fp64 : {false,true}) {
            auto code=shader("compute","cs_5_0",fp64);
            D3D12_COMPUTE_PIPELINE_STATE_DESC pd = {};
            pd.pRootSignature=root.Get(); pd.CS=bytecode(code.Get());
            dxCheck(device->CreateComputePipelineState(&pd,IID_PPV_ARGS(fp64 ? &computeD : &computeF)), "Create compute pipeline");
        }
        auto vs=shader("vertex","vs_5_0"), ps=shader("pixel","ps_5_0");
        D3D12_GRAPHICS_PIPELINE_STATE_DESC pd = {};
        pd.pRootSignature=root.Get(); pd.VS=bytecode(vs.Get()); pd.PS=bytecode(ps.Get());
        pd.BlendState.RenderTarget[0].RenderTargetWriteMask=D3D12_COLOR_WRITE_ENABLE_ALL;
        pd.RasterizerState.FillMode=D3D12_FILL_MODE_SOLID;
        pd.RasterizerState.CullMode=D3D12_CULL_MODE_NONE;
        pd.RasterizerState.DepthClipEnable=TRUE;
        pd.DepthStencilState.DepthEnable=FALSE; pd.DepthStencilState.StencilEnable=FALSE;
        pd.SampleMask=UINT_MAX; pd.PrimitiveTopologyType=D3D12_PRIMITIVE_TOPOLOGY_TYPE_TRIANGLE;
        pd.NumRenderTargets=1; pd.RTVFormats[0]=sc.Format; pd.SampleDesc.Count=1;
        dxCheck(device->CreateGraphicsPipelineState(&pd,IID_PPV_ARGS(&display)), "Create display pipeline");
        cursorReadback=buffer(32,D3D12_HEAP_TYPE_READBACK,D3D12_RESOURCE_STATE_COPY_DEST);
        resize(w,h);
        return gpu;
    }
    void resize(UINT w, UINT h) {
        if (!swap || (w==width && h==height)) return;
        wait();
        for (auto& b : back) b.Reset();
        dxCheck(swap->ResizeBuffers(2,w,h,DXGI_FORMAT_R8G8B8A8_UNORM,0), "Resize swap chain");
        auto handle=rtvHeap->GetCPUDescriptorHandleForHeapStart();
        for (UINT i=0;i<2;++i) {
            dxCheck(swap->GetBuffer(i,IID_PPV_ARGS(&back[i])), "Get back buffer");
            device->CreateRenderTargetView(back[i].Get(),nullptr,handle);
            handle.ptr+=rtvStride;
        }
        acc=buffer(UINT64(w)*h*32,D3D12_HEAP_TYPE_DEFAULT,D3D12_RESOURCE_STATE_UNORDERED_ACCESS,
                   D3D12_RESOURCE_FLAG_ALLOW_UNORDERED_ACCESS);
        width=w; height=h;
    }
    void matrices(const std::vector<float>& f, const std::vector<double>& d) {
        matF=upload(f.data(),f.size()*sizeof(float));
        matD=upload(d.data(),d.size()*sizeof(double));
    }
    void compute(const ViewConstants& view, bool fp64) {
        begin();
        commands->SetPipelineState(fp64 ? computeD.Get() : computeF.Get());
        commands->SetComputeRootSignature(root.Get());
        commands->SetComputeRoot32BitConstants(0,sizeof(view)/4,&view,0);
        commands->SetComputeRootShaderResourceView(1,(fp64 ? matD : matF)->GetGPUVirtualAddress());
        commands->SetComputeRootUnorderedAccessView(2,acc->GetGPUVirtualAddress());
        commands->Dispatch((width+15)/16,(height+15)/16,1);
        D3D12_RESOURCE_BARRIER barrier = {};
        barrier.Type=D3D12_RESOURCE_BARRIER_TYPE_UAV; barrier.UAV.pResource=acc.Get();
        commands->ResourceBarrier(1,&barrier);
        finish();
    }
    void draw(const ViewConstants& view) {
        lastBack=swap->GetCurrentBackBufferIndex();
        begin();
        transition(acc.Get(),D3D12_RESOURCE_STATE_UNORDERED_ACCESS,D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE);
        transition(back[lastBack].Get(),D3D12_RESOURCE_STATE_PRESENT,D3D12_RESOURCE_STATE_RENDER_TARGET);
        commands->SetPipelineState(display.Get());
        commands->SetGraphicsRootSignature(root.Get());
        commands->SetGraphicsRoot32BitConstants(0,sizeof(view)/4,&view,0);
        commands->SetGraphicsRootShaderResourceView(3,acc->GetGPUVirtualAddress());
        D3D12_VIEWPORT vp={0,0,(float)width,(float)height,0,1};
        D3D12_RECT rect={0,0,(LONG)width,(LONG)height};
        commands->RSSetViewports(1,&vp); commands->RSSetScissorRects(1,&rect);
        auto handle=rtvHeap->GetCPUDescriptorHandleForHeapStart(); handle.ptr+=lastBack*rtvStride;
        commands->OMSetRenderTargets(1,&handle,FALSE,nullptr);
        commands->IASetPrimitiveTopology(D3D_PRIMITIVE_TOPOLOGY_TRIANGLELIST);
        commands->DrawInstanced(3,1,0,0);
        transition(back[lastBack].Get(),D3D12_RESOURCE_STATE_RENDER_TARGET,D3D12_RESOURCE_STATE_PRESENT);
        transition(acc.Get(),D3D12_RESOURCE_STATE_PIXEL_SHADER_RESOURCE,D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
        finish();
    }
    void present(bool vsync) { dxCheck(swap->Present(vsync ? 1 : 0,0), "Present frame"); wait(); }
    void readPixel(UINT x, UINT y, double values[4]) {
        begin();
        transition(acc.Get(),D3D12_RESOURCE_STATE_UNORDERED_ACCESS,D3D12_RESOURCE_STATE_COPY_SOURCE);
        commands->CopyBufferRegion(cursorReadback.Get(),0,acc.Get(),(UINT64(y)*width+x)*32,32);
        transition(acc.Get(),D3D12_RESOURCE_STATE_COPY_SOURCE,D3D12_RESOURCE_STATE_UNORDERED_ACCESS);
        finish();
        void* data=nullptr;
        D3D12_RANGE range={0,32}, empty={0,0};
        dxCheck(cursorReadback->Map(0,&range,&data), "Read pixel");
        std::memcpy(values,data,32); cursorReadback->Unmap(0,&empty);
    }
    std::vector<uint8_t> pixels() {
        D3D12_PLACED_SUBRESOURCE_FOOTPRINT footprint;
        UINT64 bytes;
        auto desc=back[lastBack]->GetDesc();
        device->GetCopyableFootprints(&desc,0,1,0,&footprint,nullptr,nullptr,&bytes);
        auto readback=buffer(bytes,D3D12_HEAP_TYPE_READBACK,D3D12_RESOURCE_STATE_COPY_DEST);
        begin();
        transition(back[lastBack].Get(),D3D12_RESOURCE_STATE_PRESENT,D3D12_RESOURCE_STATE_COPY_SOURCE);
        D3D12_TEXTURE_COPY_LOCATION src = {}, dst = {};
        src.pResource=back[lastBack].Get(); src.Type=D3D12_TEXTURE_COPY_TYPE_SUBRESOURCE_INDEX;
        dst.pResource=readback.Get(); dst.Type=D3D12_TEXTURE_COPY_TYPE_PLACED_FOOTPRINT;
        dst.PlacedFootprint=footprint;
        commands->CopyTextureRegion(&dst,0,0,0,&src,nullptr);
        transition(back[lastBack].Get(),D3D12_RESOURCE_STATE_COPY_SOURCE,D3D12_RESOURCE_STATE_PRESENT);
        finish();
        uint8_t* data=nullptr;
        D3D12_RANGE range={0,(SIZE_T)bytes}, empty={0,0};
        dxCheck(readback->Map(0,&range,(void**)&data), "Read screenshot");
        std::vector<uint8_t> rgb(size_t(width)*height*3);
        for (UINT y=0;y<height;++y)
            for (UINT x=0;x<width;++x)
                std::memcpy(&rgb[(size_t(y)*width+x)*3],data+footprint.Offset+size_t(y)*footprint.Footprint.RowPitch+x*4,3);
        readback->Unmap(0,&empty);
        return rgb;
    }
};
