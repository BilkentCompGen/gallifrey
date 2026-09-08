#pragma once

#include <new>
#include <utility>

// NoDestructor<T> holds a T in inline (non-heap) storage; T's destructor is
// never run. For program-lifetime objects whose teardown at process exit
// would be pure overhead (e.g. a large in-memory index freed right before
// the OS reclaims the address space anyway).
//
// The T is constructed in place inside an aligned byte buffer, so the
// wrapper's own trivial destructor does not touch T.
template <class T>
class NoDestructor {
public:
    template <class... Args>
    explicit NoDestructor(Args &&...args) {
        ::new (static_cast<void *>(&space_)) T(std::forward<Args>(args)...);
    }

    NoDestructor(const NoDestructor &) = delete;
    NoDestructor &operator=(const NoDestructor &) = delete;

    T &operator*() { return *get(); }
    T *operator->() { return get(); }
    const T &operator*() const { return *get(); }
    const T *operator->() const { return get(); }

    // std::launder is required: a pointer reinterpret_cast from the byte buffer
    // is not pointer-interconvertible with the T constructed inside it.
    T *get() { return std::launder(reinterpret_cast<T *>(&space_)); }
    const T *get() const {
        return std::launder(reinterpret_cast<const T *>(&space_));
    }

private:
    alignas(T) unsigned char space_[sizeof(T)];
};
