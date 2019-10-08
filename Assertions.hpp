#ifndef ASSERTIONS_H
#define ASSERTIONS_H

#include <iostream>
#include <string>

enum class ErrType
{
 efatal,
 ewarning
};

inline static void Error(ErrType type,
			 const char *routine,
			 int lineNumber,
			 const char *msg,
			 unsigned int level)
{
    std::string baseMsg = "Level " + std::to_string(level) +
	" assertion violation\n";
    baseMsg += "Where  : " + std::string(routine) + "[" +
	std::to_string(lineNumber) + "]\n Message : ";
    baseMsg += std::string(msg);

    switch (type)
    {
    case ErrType::efatal:
	std::cout << "Fatal  : " << baseMsg << std::endl;
	std::exit(1);
	break;
    case ErrType::ewarning:
	std::cout << "Warning  : " << baseMsg << std::endl;
	break;
    }
}

#define ASSERTL0(condition, msg) \
    if (!(condition))		 \
    { \
        Error(ErrType::efatal, __FILE__, __LINE__, msg, 0);	\
    }

#endif /* ASSERTIONS_H */
