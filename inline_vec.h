#ifndef _LASM_INLINE_VEC_H_
#define _LASM_INLINE_VEC_H_

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <limits>
#include <memory>
#include <type_traits>
#include <utility>

// A vector holding up to N elements inside the object; the (N + 1)th moves
// them all into one heap block, laid out as a gtl::vector's. Graph nodes keep
// their fragments and predecessors in these, so a node with one of each owns
// no heap block. Iterators are pointers: it converts to std::span and works
// with the standard algorithms.
template <class T, uint32_t N>
class inline_vec {
   public:
    using value_type = T;
    using iterator = T *;
    using const_iterator = const T *;

    inline_vec() noexcept : heap_{nullptr} {}
    inline_vec(inline_vec &&other) noexcept : heap_{nullptr} { take(other); }
    inline_vec &operator=(inline_vec &&other) noexcept {
        if (this != &other) {
            release();
            take(other);
        }
        return *this;
    }
    inline_vec(const inline_vec &) = delete;
    inline_vec &operator=(const inline_vec &) = delete;
    ~inline_vec() { release(); }

    size_t size() const noexcept { return size_; }
    bool empty() const noexcept { return size_ == 0; }
    T *data() noexcept {
        if (spilled()) return heap_;
        return inline_;
    }
    const T *data() const noexcept {
        if (spilled()) return heap_;
        return inline_;
    }
    T *begin() noexcept { return data(); }
    T *end() noexcept { return data() + size_; }
    const T *begin() const noexcept { return data(); }
    const T *end() const noexcept { return data() + size_; }
    auto rbegin() const noexcept { return std::make_reverse_iterator(end()); }
    auto rend() const noexcept { return std::make_reverse_iterator(begin()); }
    T &operator[](size_t i) noexcept { return data()[i]; }
    const T &operator[](size_t i) const noexcept { return data()[i]; }
    T &front() noexcept { return data()[0]; }
    const T &front() const noexcept { return data()[0]; }
    T &back() noexcept { return data()[size_ - 1]; }
    const T &back() const noexcept { return data()[size_ - 1]; }

    template <class... Args>
    T &emplace_back(Args &&...args) {
        if (size_ == cap_) grow(size_t{cap_} * 2);
        T *slot =
            std::construct_at(data() + size_, std::forward<Args>(args)...);
        ++size_;
        return *slot;
    }
    void push_back(T &&v) { emplace_back(std::move(v)); }
    void push_back(const T &v) { emplace_back(v); }
    void pop_back() noexcept { std::destroy_at(data() + --size_); }
    T *erase(const T *pos) { return erase(pos, pos + 1); }
    T *erase(const T *first, const T *last) {
        T *f = begin() + (first - begin());
        T *l = begin() + (last - begin());
        T *new_end = std::move(l, end(), f);
        std::destroy(new_end, end());
        size_ -= static_cast<uint32_t>(l - f);
        return f;
    }
    // Destroys the elements and keeps the storage, as std::vector::clear.
    void clear() noexcept {
        std::destroy(begin(), end());
        size_ = 0;
    }
    void reserve(size_t n) {
        if (n > cap_) grow(n);
    }

   private:
    bool spilled() const noexcept { return cap_ > N; }
    void grow(size_t cap) {
        assert(cap <= std::numeric_limits<uint32_t>::max());
        T *block = std::allocator<T>{}.allocate(cap);
        if (spilled()) {
            std::uninitialized_move(begin(), end(), block);
        } else {
            move_inline_to(block);
        }
        std::destroy(begin(), end());
        if (spilled()) std::allocator<T>{}.deallocate(heap_, cap_);
        heap_ = block;
        cap_ = static_cast<uint32_t>(cap);
    }
    void release() noexcept {
        clear();
        if (spilled()) std::allocator<T>{}.deallocate(heap_, cap_);
        heap_ = nullptr;
        cap_ = N;
    }
    // Precondition: *this is released (inline, empty).
    void take(inline_vec &other) noexcept {
        size_ = other.size_;
        if (other.spilled()) {
            heap_ = other.heap_;
            cap_ = other.cap_;
            other.heap_ = nullptr;
            other.cap_ = N;
            other.size_ = 0;
            return;
        }
        other.move_inline_to(inline_);
        other.clear();
    }
    // Move-constructs the inline elements at dst, which holds at least N. A
    // trivially copyable T is copied as the whole buffer: GCC cannot prove
    // size_ <= N and flags the size_-element copy with -Warray-bounds.
    void move_inline_to(T *dst) noexcept {
        if constexpr (std::is_trivially_copyable_v<T>) {
            std::memcpy(dst, inline_, sizeof(inline_));
        } else {
            std::uninitialized_move(inline_, inline_ + size_, dst);
        }
    }

    uint32_t size_ = 0;
    uint32_t cap_ = N;  // N while the elements are inline
    union {
        T inline_[N];
        T *heap_;
    };
};

#endif  // _LASM_INLINE_VEC_H_
