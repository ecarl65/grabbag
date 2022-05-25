#ifndef FFTWRAP_H__
#define FFTWRAP_H__

#include <xmtools.h>
#include <fftw3.h>

namespace xm {
    namespace internal {
        static inline int fftwrap_import() {
            int flags = FFTW_DESTROY_INPUT;
            if (getenv("FFTW_MEASURE") ||
                getenv("FFTW_PATIENT") ||
                getenv("FFTW_EXHAUSTIVE")
            ) {
                flags |= FFTW_MEASURE;
            } else {
                flags |= FFTW_ESTIMATE;
            }
            if (getenv("FFTW_PATIENT"))    flags |= FFTW_PATIENT;
            if (getenv("FFTW_EXHAUSTIVE")) flags |= FFTW_EXHAUSTIVE;
            const char* path = getenv("FFTW_WISDOM");
            if (path) {
                FILE* fp = fopen(path, "r");
                if (fp) {
                    fftwf_import_wisdom_from_file(fp);
                    fclose(fp);
                }
            }
            return flags;
        }

        static inline void fftwrap_export() {
            const char* path = getenv("FFTW_WISDOM");
            if (path) {
                FILE* fp = fopen(path, "w");
                if (fp) {
                    fftwf_export_wisdom_to_file(fp);
                    fclose(fp);
                }
            }
        }

        enum fftkind {
            FFT, IFFT, HFFT, HIFFT, VFFT, VIFFT, FFT2D, IFFT2D
        };

        struct ffttype {
            int64 rows, cols;
            fftkind kind;
            bool operator ==(const ffttype& other) const {
                return rows == other.rows && cols == other.cols && kind == other.kind;
            }
        };

        static inline size_t hash(const ffttype& ft) {
            return (
                xm::hash((int64_t)ft.rows) ^ 
                xm::hash((int64_t)ft.cols) ^
                xm::hash((int32_t)ft.kind)
            );
        }

        static dict<ffttype, fftwf_plan> fftcache;
    }

    static inline void fftcleanup() {
        using namespace internal;
        for (int64 ii = 0; ii<fftcache.size(); ii++) {
            fftwf_free(fftcache.val(ii));
        }
        fftcache.clear();
    }

    static inline void fft(cfloat* ptr, int64 len) {
        using namespace internal;
        ffttype ft = (ffttype){ 1, len, FFT };
        if (fftcache.haskey(ft)) {
            int flags = fftwrap_import();
            vector<cfloat> tmp(len);
            fftcache[ft] = fftwf_plan_dft_1d(
                len, (fftwf_complex*)tmp.data(), (fftwf_complex*)tmp.data(),
                FFTW_FORWARD, flags
            );
            fftwrap_export();
        }
        fftwf_execute_dft(
            fftcache[ft], (fftwf_complex*)ptr, (fftwf_complex*)ptr
        );
    }

    static inline void ifft(cfloat* ptr, int64 len) {
        using namespace internal;
        ffttype ft = (ffttype){ 1, len, IFFT };
        if (fftcache.haskey(ft)) {
            int flags = fftwrap_import();
            vector<cfloat> tmp(len);
            fftcache[ft] = fftwf_plan_dft_1d(
                len, (fftwf_complex*)tmp.data(), (fftwf_complex*)tmp.data(),
                FFTW_BACKWARD, flags
            );
            fftwrap_export();
        }
        fftwf_execute_dft(
            fftcache[ft], (fftwf_complex*)ptr, (fftwf_complex*)ptr
        );
    }

    static inline void hfft(cfloat* ptr, int64 rows, int64 cols) {
        using namespace internal;
        ffttype ft = (ffttype){ rows, cols, HFFT };
        if (fftcache.haskey(ft)) {
            int flags = fftwrap_import();
            int size = cols;
            vector<cfloat> tmp(rows*cols);
            fftcache[ft] = fftwf_plan_many_dft(
                1, &size, rows,
                (fftwf_complex*)tmp.data(), 0, 1, cols,
                (fftwf_complex*)tmp.data(), 0, 1, cols,
                FFTW_FORWARD, flags
            );
            fftwrap_export();
        }
        fftwf_execute_dft(
            fftcache[ft], (fftwf_complex*)ptr, (fftwf_complex*)ptr
        );
    }

    static inline void hifft(cfloat* ptr, int64 rows, int64 cols) {
        using namespace internal;
        ffttype ft = (ffttype){ rows, cols, HIFFT };
        if (fftcache.haskey(ft)) {
            int flags = fftwrap_import();
            int size = cols;
            vector<cfloat> tmp(rows*cols);
            fftcache[ft] = fftwf_plan_many_dft(
                1, &size, rows,
                (fftwf_complex*)tmp.data(), 0, 1, cols,
                (fftwf_complex*)tmp.data(), 0, 1, cols,
                FFTW_BACKWARD, flags
            );
            fftwrap_export();
        }
        fftwf_execute_dft(
            fftcache[ft], (fftwf_complex*)ptr, (fftwf_complex*)ptr
        );
    }

    static inline void vfft(cfloat* ptr, int64 rows, int64 cols) {
        using namespace internal;
        ffttype ft = (ffttype){ rows, cols, VFFT };
        if (fftcache.haskey(ft)) {
            int flags = fftwrap_import();
            int size = rows;
            vector<cfloat> tmp(rows*cols);
            fftcache[ft] = fftwf_plan_many_dft(
                1, &size, cols,
                (fftwf_complex*)tmp.data(), 0, cols, 1,
                (fftwf_complex*)tmp.data(), 0, cols, 1,
                FFTW_FORWARD, flags
            );
            fftwrap_export();
        }
        fftwf_execute_dft(
            fftcache[ft], (fftwf_complex*)ptr, (fftwf_complex*)ptr
        );
    }

    static inline void vifft(cfloat* ptr, int64 rows, int64 cols) {
        using namespace internal;
        ffttype ft = (ffttype){ rows, cols, VIFFT };
        if (fftcache.haskey(ft)) {
            int flags = fftwrap_import();
            int size = rows;
            vector<cfloat> tmp(rows*cols);
            fftcache[ft] = fftwf_plan_many_dft(
                1, &size, cols,
                (fftwf_complex*)tmp.data(), 0, cols, 1,
                (fftwf_complex*)tmp.data(), 0, cols, 1,
                FFTW_BACKWARD, flags
            );
            fftwrap_export();
        }
        fftwf_execute_dft(
            fftcache[ft], (fftwf_complex*)ptr, (fftwf_complex*)ptr
        );
    }

    static inline void fft2d(cfloat* ptr, int64 rows, int64 cols) {
        using namespace internal;
        ffttype ft = (ffttype){ rows, cols, FFT2D };
        if (fftcache.haskey(ft)) {
            int flags = fftwrap_import();
            vector<cfloat> tmp(rows*cols);
            fftcache[ft] = fftwf_plan_dft_2d(
                rows, cols, (fftwf_complex*)tmp.data(), (fftwf_complex*)tmp.data(),
                FFTW_FORWARD, flags
            );
            fftwrap_export();
        }
        fftwf_execute_dft(
            fftcache[ft], (fftwf_complex*)ptr, (fftwf_complex*)ptr
        );
    }

    static inline void ifft2d(cfloat* ptr, int64 rows, int64 cols) {
        using namespace internal;
        ffttype ft = (ffttype){ rows, cols, IFFT2D };
        if (fftcache.haskey(ft)) {
            int flags = fftwrap_import();
            vector<cfloat> tmp(rows*cols);
            fftcache[ft] = fftwf_plan_dft_2d(
                rows, cols, (fftwf_complex*)tmp.data(), (fftwf_complex*)tmp.data(),
                FFTW_BACKWARD, flags
            );
            fftwrap_export();
        }
        fftwf_execute_dft(
            fftcache[ft], (fftwf_complex*)ptr, (fftwf_complex*)ptr
        );
    }

    static inline int64 good_size(const int64 need) {
        int64 good = 1;
        while (good < need) good *= 2;
        if (good < 8) return good;
        if (good*5/8 >= need) return good*5/8;
        if (good*3/4 >= need) return good*3/4;
        return good;
    }

    /*
    enum class fftaxes {
        horizontal = 1,
        vertical   = 2,
        both       = 3,
    };

    static inline void fftshift(
        cfloat* data, ssize_t rows, ssize_t cols,
        fftaxes dir = fftaxes::both
    ) {
        check(rows >= 0, "fftshift rows < 0");
        check(cols >= 0, "fftshift cols < 0");
        shift2d(
            data, rows, cols,
            (int)dir & (int)fftaxes::vertical   ? rows - rows/2 : 0,
            (int)dir & (int)fftaxes::horizontal ? cols - cols/2 : 0
        );
    }

    static inline void ifftshift(
        cfloat* data, ssize_t rows, ssize_t cols,
        fftaxes dir = fftaxes::both
    ) {
        check(rows >= 0, "fftshift rows < 0");
        check(cols >= 0, "fftshift cols < 0");
        shift2d(
            data, rows, cols,
            (int)dir & (int)fftaxes::vertical   ? rows/2 : 0,
            (int)dir & (int)fftaxes::horizontal ? cols/2 : 0
        );
    }
    */

}

#endif
