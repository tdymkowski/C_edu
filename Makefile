CC := gcc
CFLAGS :=-Wall -Wextra -std=c11 -g
LDFLAGS := -lm

SRCS := main.c xyz_reader.c lj_potential.c

HEADERS := xyz_reader.h lj_potential.h main.h

OBJS=$(SRCS:.c=.o)

TARGET=lj_potential

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.c $(HEADERS)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)
