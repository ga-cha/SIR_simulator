function [] = test_write(i)
    write_async(array2table(ones(10000,20)*i), 'test.csv')
end