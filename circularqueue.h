#pragma once

#ifndef CIRCULAR_QUEUE_H
#define CIRCULAR_QUEUE_H
#include <stddef.h>
#include <optional>
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


template<typename T, size_t Capacity>
class cdqa_t {
    static_assert(Capacity > 0, "Capacity must be greater than 0");

private:
    T buffer[Capacity];
    size_t head = 0;
    size_t tail = 0;
    size_t count = 0;

    size_t wrap(size_t index) const {
        return index % Capacity;
    }

public:
    cdqa_t() = default;

    size_t size() const { return count; }
    constexpr size_t capacity() const { return Capacity; }
    bool empty() const { return count == 0; }
    bool full() const { return count == Capacity; }

    bool push_front(const T& value) {
        if (full()) return false;
        head = (head + Capacity - 1) % Capacity;
        buffer[head] = value;
        ++count;
        return true;
    }

    bool push_back(const T& value) {
        if (full()) return false;
        buffer[tail] = value;
        tail = (tail + 1) % Capacity;
        ++count;
        return true;
    }

    bool pop_front() {
        if (empty()) return false;
        head = (head + 1) % Capacity;
        --count;
        return true;
    }

    bool pop_back() {
        if (empty()) return false;
        tail = (tail + Capacity - 1) % Capacity;
        --count;
        return true;
    }

    std::optional<T> front() const {
        if (empty()) return std::nullopt;
        return buffer[head];
    }

    std::optional<T> back() const {
        if (empty()) return std::nullopt;
        return buffer[(tail + Capacity - 1) % Capacity];
    }

    std::optional<T> at(size_t index) const {
        if (index >= count) return std::nullopt;
        return buffer[wrap(head + index)];
    }

    T& operator[](size_t index) {
        return buffer[wrap(head + index)];
    }

    const T& operator[](size_t index) const {
        return buffer[wrap(head + index)];
    }

    void clear() {
        head = 0;
        tail = 0;
        count = 0;
    }
};


template<class T, int N>
struct smq{
    cdqa_t<T, N> miniq;
    cdqa_t<T, N> fifoq;
    void push_back(T elem){
        fifoq.push_back(elem);
        if(miniq.empty() == false && elem < miniq.front()){
            miniq.clear();
        }
        while(miniq.empty() == false && elem < miniq.back()){
            miniq.pop_back();
        }
        miniq.push_back(elem);
    }   
    void pop_front(){
        T elem = fifoq.front().value();
        fifoq.pop_front();

        if (elem == miniq.front()){
            miniq.pop_front();
        }
    }
    constexpr T front() const{
        return fifoq.front().value();
    }
    constexpr T back() const {
        return fifoq.back().value();
    }
    constexpr T min() const {
        return miniq.front().value();
    }
    T pushpop(T elem){
        push_back(elem);
        elem = front();
        pop_front();
        return elem;
    }
    T& operator[](std::size_t i) {
        return fifoq[i].value();
    }

    constexpr const T& operator[](std::size_t i) const {
        return fifoq[i].value();
    }

    constexpr std::size_t size() const {
        return fifoq.size();
    }

    constexpr bool empty() const {
        return fifoq.empty();
    }
};

#endif
