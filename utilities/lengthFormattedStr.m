function str = lengthFormattedStr(str, maxChar)
% Properly format the size of a string before writing to file
        lenChar = length(str); 
        str = [str(1:min(lenChar,maxChar)) blanks(maxChar-lenChar + 6)];
end