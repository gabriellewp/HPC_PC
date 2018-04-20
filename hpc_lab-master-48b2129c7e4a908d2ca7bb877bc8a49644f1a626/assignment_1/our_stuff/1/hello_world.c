#include <stdio.h>
#include <unistd.h>

int main(int argc, char** argv)
{
	char hostname[1024];
	gethostname(hostname, 1024);
	printf("Hello World! (from: %s)\n", hostname);

	return 0;
}
