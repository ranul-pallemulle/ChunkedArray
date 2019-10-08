#ifndef ARRAYPOLICIES_HPP
#define ARRAYPOLICIES_HPP

template<typename ObjectType, typename enabled = void>
class ArrayInitializationPolicy;

template<typename ObjectType>
class ArrayInitializationPolicy<ObjectType> {
    
public:
    
    static void Initialize(ObjectType* data, unsigned int itemsToCreate)
    {
	DoInitialization(
			 data, itemsToCreate,
			 [](ObjectType *element) { new (element) ObjectType;}
			 );
    }

    static void Initialize(ObjectType* data, unsigned int itemsToCreate, const ObjectType& initValue)
    {
	DoInitialization(
			 data, itemsToCreate,
			 [&](ObjectType *element) { new (element) ObjectType(initValue); 
			 });
    }


private:
    template<typename CreateType>
    static void DoInitialization(ObjectType* data, unsigned int itemsToCreate, const CreateType &f)
    {
	unsigned int nextObjectToCreate = 0;
	try {
	    for (unsigned int i = 0; i < itemsToCreate; ++i) {
		ObjectType* memLocation = &data[i];
		f(memLocation);
		++nextObjectToCreate;
	    }
	}
	catch(...) {
	    for (unsigned int i = 0; i < nextObjectToCreate; ++i) {
		ObjectType* memLocation = &data[nextObjectToCreate - i - 1];
		memLocation->~ObjectType();
	    }
	    throw;
	}
    }
};

template<typename ObjectType, typename enabled = void>
class ArrayDestructionPolicy;

template<typename ObjectType>
class ArrayDestructionPolicy<ObjectType> {
public:
    static void Destroy(ObjectType* data, unsigned int itemsToDestroy)
    {
	for (unsigned int i = 0; i < itemsToDestroy; ++i) {
	    ObjectType* memLocation = &data[itemsToDestroy - i - 1];
	    memLocation->~ObjectType();
	}
    }
};

#endif
