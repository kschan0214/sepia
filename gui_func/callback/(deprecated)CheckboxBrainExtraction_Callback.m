function CheckboxBrainExtraction_Callback(source,eventdata,h)
% if BET checkbox is checked then empty mask edit field and disable open
% pushbutton

% global h

if ~h.dataIO.checkbox.brainExtraction.Value
    set(h.dataIO.button.maskdir,        'Enable','on');
    set(h.dataIO.edit.maskdir,          'Enable','on');
    set(h.dataIO.popup.brainExtraction, 'Enable', 'off');
    set(h.dataIO.edit.fractionalThres,  'Enable', 'off');
    set(h.dataIO.edit.gradientThres,    'Enable', 'off')
else
    set(h.dataIO.button.maskdir,        'Enable','off');
    set(h.dataIO.edit.maskdir,          'Enable','off');
    set(h.dataIO.edit.maskdir,          'String','');
    set(h.dataIO.popup.brainExtraction, 'Enable', 'on');
    set(h.dataIO.edit.fractionalThres,  'Enable', 'on');
    set(h.dataIO.edit.gradientThres,    'Enable', 'on');
end

end