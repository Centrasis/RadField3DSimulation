#include <vector>

template <typename T>
class RollingBuffer {
public:
    RollingBuffer(size_t maxSize) : pointer(0), buffer(maxSize) {}

    void push(const T& value) {
        if (!this->is_full && pointer + 1 == this->buffer.size())
            this->is_full = true;
        buffer[pointer] = value;
        pointer = ++pointer % this->buffer.size();
    }

    const T& operator[](size_t index) const {
        return buffer[index];
    }

    size_t size() const {
        return (is_full) ? buffer.size() : pointer;
    }

    void reset() {
        pointer = 0;
        is_full = false;
    }

private:
    std::vector<T> buffer;
    size_t pointer;
    bool is_full = false;
};