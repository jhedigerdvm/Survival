<?xml version="1.0" encoding="UTF-8"?>
<ArcPad>
   <LAYER name="Reservoir">
      <SYMBOLOGY>
         <SIMPLELABELRENDERER visible="false" field="WELL_NAME">
            <TEXTSYMBOL angle="0" fontcolor="0,0,0" font="Arial" fontsize="8" horzalignment="left" vertalignment="bottom" rtl="false">
            </TEXTSYMBOL>
         </SIMPLELABELRENDERER>
         <VALUEMAPRENDERER lookupfield="TYPE">
            <EXACT value=" " label=" ">
               <SIMPLEMARKERSYMBOL angle="0" width="4" color="139,166,133" type="circle" outline="0,0,0" outlinewidth="1" outlinetype="solid">
               </SIMPLEMARKERSYMBOL>
            </EXACT>
            <EXACT value="Pila" label="Pila">
               <TRUETYPEMARKERSYMBOL font="ESRI Default Marker" character="46" fontcolor="0,112,255" fontsize="18" angle="0">
               </TRUETYPEMARKERSYMBOL>
            </EXACT>
            <EXACT value="Trough" label="Trough">
               <GROUPSYMBOL>
                  <TRUETYPEMARKERSYMBOL font="ESRI Default Marker" character="46" fontcolor="0,0,0" fontsize="12" angle="0">
                  </TRUETYPEMARKERSYMBOL>
                  <TRUETYPEMARKERSYMBOL font="ESRI Default Marker" character="33" fontcolor="0,112,255" fontsize="12" angle="0">
                  </TRUETYPEMARKERSYMBOL>
               </GROUPSYMBOL>
            </EXACT>
            <OTHER label="&lt;all other values&gt;">
               <SIMPLEMARKERSYMBOL angle="0" width="4" color="0,143,10" type="circle" outline="0,0,0" outlinewidth="1" outlinetype="solid">
               </SIMPLEMARKERSYMBOL>
            </OTHER>
         </VALUEMAPRENDERER>
      </SYMBOLOGY>
   </LAYER>
</ArcPad>
