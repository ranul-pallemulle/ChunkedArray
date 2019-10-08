#ifndef MEMORYMANAGER_HPP
#define MEMORYMANAGER_HPP

template<typename DataType>
class MemoryManager
{
public:
    static DataType* RawAllocate(unsigned int NumberOfElements)
    {
	return static_cast<DataType*>(::operator new(NumberOfElements * sizeof(DataType)));
    }

    static void RawDeallocate(DataType* array, unsigned int)
    {
	::operator delete(array);
    }
	
};

#endif //MEMORYMANAGER_HPP
