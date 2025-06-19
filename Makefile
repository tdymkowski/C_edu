CC := gcc
CFLAGS :=-Wall -Wextra -std=c11 -g
LDFLAGS := -lm

SRCS := main.c xyz_reader.c xyz_writer.c lj_potential.c read_parameters.c verlet.c

HEADERS := read_parameters.h xyz_reader.h xyz_writer.h lj_potential.h verlet.h main.h 

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
