.PHONY: all test check-release clean

all:
	$(MAKE) -C src

test:
	$(MAKE) -C src test
	sh scripts/check_release.sh

check-release:
	sh scripts/check_release.sh

clean:
	$(MAKE) -C src clean
