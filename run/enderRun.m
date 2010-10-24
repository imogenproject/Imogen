function enderRun()

%=== ENDERRUN ======================================================================================
% Test if the program is running on linux/unix. If so, end the Matlab process or if not return
% operation complete value.
%===================================================================================================

    if ( isunix() == true )    
        exit
    else
        disp( 'Operation Complete' );
    end
    
end