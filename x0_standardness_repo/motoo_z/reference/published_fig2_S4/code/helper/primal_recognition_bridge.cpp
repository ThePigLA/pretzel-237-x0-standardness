// Recognition-only bridge for the exact pinned Regina 7.4.1 Linux wheel.
// No Regina object field is inspected or modified. The only engine calls
// are its exported public tightEncode(), and original isSphere()/isBall().
#include <atomic>
#include <cstddef>
#include <dlfcn.h>
#include <sstream>
#include <stdexcept>
#include <string>

using Recognizer = bool (*)(const void*);
using Encoder = void (*)(const void*, std::ostream&);
using ProofCallback = int (*)(const char*, std::size_t, int);
static Recognizer original_sphere = nullptr;
static Recognizer original_ball = nullptr;
static Encoder encode = nullptr;
static ProofCallback proof_callback = nullptr;
static void* engine_handle = nullptr;
static std::string configuration_error;
static std::atomic<unsigned long long> counts[5]{};
static thread_local unsigned depth = 0;
static const char* sphere_symbol = "_ZNK6regina13TriangulationILi3EE8isSphereEv";
static const char* ball_symbol = "_ZNK6regina13TriangulationILi3EE6isBallEv";
static const char* encoder_symbol = "_ZNK6regina6detail17TriangulationBaseILi3EE11tightEncodeERSo";

extern "C" bool sphere_bridge(const void*)
    asm("_ZNK6regina13TriangulationILi3EE8isSphereEv");
extern "C" bool ball_bridge(const void*)
    asm("_ZNK6regina13TriangulationILi3EE6isBallEv");

extern "C" int motoo_primal_configure(const char* path, ProofCallback callback) {
    try {
        if (!path || !callback)
            throw std::runtime_error("missing engine path or proof callback");
        void* handle = dlopen(path, RTLD_NOW | RTLD_NOLOAD);
        if (!handle)
            throw std::runtime_error(dlerror());
        auto sphere = reinterpret_cast<Recognizer>(dlsym(handle, sphere_symbol));
        auto ball = reinterpret_cast<Recognizer>(dlsym(handle, ball_symbol));
        auto serializer = reinterpret_cast<Encoder>(dlsym(handle, encoder_symbol));
        if (!sphere || !ball || !serializer)
            throw std::runtime_error("required public engine symbol is missing");
        if (sphere == sphere_bridge || ball == ball_bridge)
            throw std::runtime_error("original recognition symbol resolved to bridge");
        engine_handle = handle;
        original_sphere = sphere;
        original_ball = ball;
        encode = serializer;
        proof_callback = callback;
        return 1;
    } catch (const std::exception& e) {
        configuration_error = e.what();
        return 0;
    } catch (...) {
        configuration_error = "unknown configuration error";
        return 0;
    }
}

extern "C" const char* motoo_primal_error() {
    return configuration_error.c_str();
}

extern "C" unsigned long long motoo_primal_counter(unsigned index) {
    return index < 5 ? counts[index].load() : 0;
}

static bool recognize(const void* object, int kind) {
    if (!object)
        throw std::invalid_argument("null triangulation passed to recognition");
    counts[kind].fetch_add(1);
    if (proof_callback && encode && depth == 0) {
        ++depth;
        int verdict = 0;
        try {
            // Triangulation<3> has a single non-virtual TriangulationBase<3>
            // base at offset zero in the pinned Itanium ABI. tightEncode()
            // is invoked as an exported method, with no field/layout access.
            std::ostringstream stream;
            encode(object, stream);
            const std::string bytes = stream.str();
            verdict = proof_callback(bytes.data(), bytes.size(), kind);
        } catch (...) {
            counts[4].fetch_add(1);
        }
        --depth;
        if (verdict == 1) {
            counts[2].fetch_add(1);
            return true;
        }
        if (verdict < 0 || verdict > 1)
            counts[4].fetch_add(1);
    }
    counts[3].fetch_add(1);
    Recognizer original = kind == 0 ? original_sphere : original_ball;
    if (!original) {
        original = reinterpret_cast<Recognizer>(
            dlsym(RTLD_NEXT, kind == 0 ? sphere_symbol : ball_symbol));
    }
    if (!original || original == sphere_bridge || original == ball_bridge)
        throw std::runtime_error("recognition bridge has no original fallback");
    // Unknown means the exact original routine decides. It never means false.
    return original(object);
}

extern "C" bool sphere_bridge(const void* object) { return recognize(object, 0); }
extern "C" bool ball_bridge(const void* object) { return recognize(object, 1); }
