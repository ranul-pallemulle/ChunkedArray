#ifndef SHARED_ARRAY_HPP
#define SHARED_ARRAY_HPP

#include "ArrayPolicies.hpp"
#include "MemoryManager.hpp"

struct OneD
{
    static const unsigned int Value = 1;
};

template<typename Dim, typename DataType>
class Array;

template<typename DataType>
class Array<OneD, const DataType>
{
public:
    Array():
	m_size(0),
	m_capacity(0),
	m_data(nullptr),
	m_count(nullptr),
	m_offset(0)
    {
	CreateStorage(m_capacity);
    }

    explicit Array(unsigned int dim1Size):
	m_size(dim1Size),
	m_capacity(dim1Size),
	m_data(nullptr),
	m_count(nullptr),
	m_offset(0)
    {
	CreateStorage(m_capacity);
	ArrayInitializationPolicy<DataType>::Initialize(m_data, m_capacity);
    }

    Array(unsigned int dim1Size, const DataType& initValue) :
	m_size(dim1Size),
	m_capacity(dim1Size),
	m_data(nullptr),
	m_count(nullptr),
	m_offset(0)
    {
	CreateStorage(m_capacity);
	ArrayInitializationPolicy<DataType>::Initialize(m_data, m_capacity, initValue);
    }

    ~Array()
    {
	if (m_count == nullptr)
	    return;
	*m_count -= 1;
	if (*m_count == 0) {
	    ArrayDestructionPolicy<DataType>::Destroy(m_data, m_capacity);
	    MemoryManager<DataType>::RawDeallocate(m_data, m_capacity);
	    delete m_count;
	}
    }

    Array<OneD, const DataType>& operator=(const Array<OneD, const DataType>& rhs)
    {
	*m_count -= 1;
	if (*m_count == 0) {
	    ArrayDestructionPolicy<DataType>::Destroy(m_data, m_capacity);
	    MemoryManager<DataType>::RawDeallocate(m_data, m_capacity);
	    delete m_count;
	}
	m_data = rhs.m_data;
	m_capacity = rhs.m_capacity;
	m_count = rhs.m_count;
	*m_count += 1;
	m_offset = rhs.m_offset;
	m_size = rhs.m_size;
	return *this;
    }

    const DataType& operator[](unsigned int i) const
    {
	return *(m_data + i + m_offset);
    }

    const DataType* get() const
    {
	return m_data + m_offset;
    }

    unsigned int num_elements() const
    {
	return m_size;
    }

protected:
    unsigned int m_size;
    unsigned int m_capacity;
    DataType* m_data;
    unsigned int* m_count;
    unsigned int m_offset;

private:
    void CreateStorage(unsigned int size)
    {
	DataType* storage = MemoryManager<DataType>::RawAllocate(size);
	m_data = storage;
	m_count = new unsigned int();
	*m_count = 1;
    }

};

template<typename DataType>
class Array<OneD, DataType> : public Array<OneD, const DataType> {
public:
    typedef Array<OneD, const DataType> BaseType;

public:
    Array() : BaseType() {}
    
    explicit Array(unsigned int dim1Size) : BaseType(dim1Size) {}

    Array(unsigned int dim1Size, const DataType& initValue) :
	BaseType(dim1Size, initValue)
    {
    }
    
    Array<OneD, DataType>& operator=(const Array<OneD, DataType>& rhs)
    {
	BaseType::operator=(rhs);
	return *this;
    }
    
    using BaseType::operator[];
    DataType& operator[](unsigned int i)
    {
	return (get())[i];
    }

    using BaseType::get;
    DataType* get(){ return this->m_data + this->m_offset;}
};


#endif //SHARED_ARRAY_HPP
