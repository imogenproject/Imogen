function partneredCloseCB(src, event, hPartner)

    if ishandle(hPartner); delete(hPartner); end
    delete(src);
    
end