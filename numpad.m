function num = numpad(num,len)
    if isnumeric(num)
        num = num2str(num);
    end
    
    if length(num) < len
        num = numpad(['0',num],len);
    end
end