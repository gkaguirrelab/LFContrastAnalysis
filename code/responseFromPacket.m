[responses,plot] =  responseFromPacket(obj, params, packetPocket)




for ii = 1:length(packetPocket)
    
    
    response.values   = obj.computeResponse(params,packetPocket{ii}.stimulus,packetPocket{ii}.kernel);
    response.timebase = packetPocket{ii}.stimulus.timebase;
    
    if createPlotStruct
        plotstruct.response = 
        plotStruct.timebase =
        plotStruct.color    = 
        
      
    
end