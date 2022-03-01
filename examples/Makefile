APP_NAME = sqrt3
OBJS += sqrt3.o
OBJS += main.o

CXX = mpicc
CXXFLAGS = -I. -O3 -Wall -Wextra

default: $(APP_NAME)

$(APP_NAME): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)

%.o: %.c
	$(CXX) $< $(CXXFLAGS) -c -o $@

clean:
	/bin/rm -rf *~ *.o $(APP_NAME)
