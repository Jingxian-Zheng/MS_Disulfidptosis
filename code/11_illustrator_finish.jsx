#target illustrator
(function(){
var root=$.getenv('REVISION_ROOT');
if(!root){var folder=Folder.selectDialog('Select the prepared reproduction workspace');if(!folder){throw new Error('No workspace selected');}root=folder.fsName;}
root=root.split(String.fromCharCode(92)).join('/');
var names=['Fig1','Fig2','Fig3','Fig4','Fig5','Fig6','S1_Fig','S2_Fig','S3_Fig','S4_Fig'];
var log=new File(root+'/figures/illustrator/Illustrator_production_log.tsv');log.encoding='UTF-8';log.open('w');log.writeln('file\tIllustrator_version\twidth_pt\theight_pt\ttext_frames\tfont_edits\tstatus');
var old=app.userInteractionLevel;app.userInteractionLevel=UserInteractionLevel.DONTDISPLAYALERTS;
try{for(var i=0;i<names.length;i++){
 var n=names[i],doc=null;
 try{
  doc=app.open(new File(root+'/figures/source/'+n+'.pdf'));
  doc.layers[0].name='Data graphics and editable text';
  var edits=0,nt=doc.textFrames.length;
  for(var j=0;j<nt;j++){
   var t=doc.textFrames[j],a=t.textRange.characterAttributes,oldname=a.textFont.name;
   var bold=/bold/i.test(oldname);
   a.textFont=app.textFonts.getByName(bold?'Arial-BoldMT':'ArialMT');
   if(a.size<8)a.size=8;
   edits++;
  }
  var nr=doc.artboards[0].artboardRect;
  var ai=new IllustratorSaveOptions();ai.pdfCompatible=true;ai.compressed=true;
  doc.saveAs(new File(root+'/figures/illustrator/'+n+'.ai'),ai);
  var png=new ExportOptionsPNG24();png.antiAliasing=true;png.transparency=false;png.artBoardClipping=true;png.horizontalScale=416.6666667;png.verticalScale=416.6666667;
  doc.exportFile(new File(root+'/figures/illustrator/'+n),ExportType.PNG24,png);
  var tif=new ExportOptionsTIFF();tif.resolution=(n=='S1_Fig'?450:600);tif.imageColorSpace=ImageColorSpace.RGB;tif.lZWCompression=true;tif.antiAliasing=AntiAliasingMethod.ARTOPTIMIZED;tif.artBoardClipping=true;
  doc.exportFile(new File(root+'/figures/illustrator/'+n),ExportType.TIFF,tif);
  if(i<6){var eps=new EPSSaveOptions();eps.embedLinkedFiles=true;eps.includeDocumentThumbnails=true;doc.saveAs(new File(root+'/figures/illustrator/'+n+'.eps'),eps);}
  var pdf=new PDFSaveOptions();pdf.preserveEditability=false;pdf.viewAfterSaving=false;
  doc.saveAs(new File(root+'/figures/illustrator/'+n+'.pdf'),pdf);
  log.writeln(n+'\t'+app.version+'\t'+(nr[2]-nr[0])+'\t'+(nr[1]-nr[3])+'\t'+nt+'\t'+edits+'\tOK');
 }catch(e){log.writeln(n+'\t'+app.version+'\tERROR\t'+e.toString());}
 finally{if(doc)doc.close(SaveOptions.DONOTSAVECHANGES);}
}}finally{app.userInteractionLevel=old;log.close();}
return 'Illustrator figure production finished';
})();
