test_dir = "/home/karol/Documents/UNIVERSITA/POLITO/shimmer/shimmer++/build/unit_tests";
disma_pin = load(strcat(test_dir , "/var_time_DISMA_pin.dat"))';
denerg_pin = load(strcat(test_dir , "/var_time_DENERG_pin.dat"))';

disp(" * Relative Error")
(abs(denerg_pin - disma_pin)./Xdenerg_pin);

