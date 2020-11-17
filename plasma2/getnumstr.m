function str = getnumstr(num,length)

str = '';
if num == 0
    for ii = 1:length
        str = [str '0'];
    end
    return
else
    for ii = 0:length
        if num < 10^ii
            for jj = 1:length-ii
                str = [str '0'];
            end
            str = [str num2str(num)];
            return
        end
    end
end