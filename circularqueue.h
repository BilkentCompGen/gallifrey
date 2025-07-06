#pragma once

#ifndef CIRCULAR_QUEUE_H
#define CIRCULAR_QUEUE_H
#include <stddef.h>
template <typename T, size_t N>
class CircularQueue {
private:
    T data[N];
    size_t frontIndex = 0;
    size_t rearIndex = N - 1;
    size_t count = 0;

public:
    constexpr bool isFull() const noexcept {
        return count == N;
    }

    constexpr bool isEmpty() const noexcept {
        return count == 0;
    }

    constexpr void enqueue(const T& value) noexcept {
        if (isFull()) {
            frontIndex = (frontIndex + 1) % N;
            --count;
        }
        rearIndex = (rearIndex + 1) % N;
        data[rearIndex] = value;
        ++count;
    }

    constexpr void dequeue() noexcept {
        if (!isEmpty()) {
            frontIndex = (frontIndex + 1) % N;
            --count;
        }
    }

    constexpr T& front() noexcept {
        return data[frontIndex];
    }

    constexpr const T& front() const noexcept {
        return data[frontIndex];
    }

    constexpr T& rear() noexcept {
        return data[rearIndex];
    }

    constexpr const T& rear() const noexcept {
        return data[rearIndex];
    }

    constexpr T& operator[](size_t idx) noexcept {
        return data[(frontIndex + idx) % N];
    }

    constexpr const T& operator[](size_t idx) const noexcept {
        return data[(frontIndex + idx) % N];
    }

    constexpr size_t size() const noexcept {
        return count;
    }

    constexpr size_t max_size() const noexcept {
        return N;
    }
};
#endif
