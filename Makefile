

C_FILES := $(wildcard Chipmunk-7.0.1/src/*.c)
OBJ_DIR := obj_c99
OBJ_FILES := $(addprefix $(OBJ_DIR)/,$(notdir $(C_FILES:.c=.o)))
LD_FLAGS := ...
CC_FLAGS := -c -DNDEBUG -IChipmunk-7.0.1/include -Wall -O3 -std=c99 -pthread
# CC_FLAGS := -c  -IChipmunk-7.0.1/include -Wall -O -std=c99 -pthread 


all: $(OBJ_FILES)

$(OBJ_DIR)/%.o: Chipmunk-7.0.1/src/%.c
#	export CXX=/usr/local/Cellar/gcc/7.3.0_1/bin/gcc-7
	gcc $(CC_FLAGS) -o $@ $<

archive:	
	mkdir -p lib
	ar rcs ./lib/libcp.a ./$(OBJ_DIR)/*.o 
clean: 
	rm -f $(OBJ_DIR)/*.o
	rm -f ./lib/libcp.a
