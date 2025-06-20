CC := gcc
CFLAGS :=-Wall -Wextra -std=c11 -g
LDFLAGS := -lm

SRCS_DIRS := . models neighbors readers writers

SRCS := $(wildcard $(addsuffix /*.c, $(SRCS_DIRS)))
HEADERS := $(wildcard $(addsuffix /*.h, $(SRCS_DIRS)))

OBJS := $(patsubst %.c, %.o, $(SRCS))

TARGET=simple_md

.PHONY: all clean

INCLUDES := $(addprefix -I, $(SRCS_DIRS))
CFLAGS += $(INCLUDES)


all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
