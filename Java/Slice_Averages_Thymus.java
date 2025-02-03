/** Computes slice-wise averages of non-zero pixel values (usually masked by an ROI) 
 over a stack of images. 
 For use by 3T NIA.
 @author: Sokratis Makrogiannis
*/

import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.filter.*;

public class Slice_Averages_Thymus implements PlugInFilter {

 protected ImageStack stack;

 public int setup(String arg, ImagePlus imp) {
  stack = imp.getStack();
  return DOES_32+STACK_REQUIRED;
 }

 public void run(ImageProcessor ip) {
  float[] pixels;
  int dimension = stack.getWidth()*stack.getHeight();
  int nSlices = stack.getSize();
  float[] sum = new float[nSlices];
  float[] count= new float[nSlices ];
  float[] average= new float[nSlices];

  for (int i=1;i<=stack.getSize();i++) {
   pixels = (float[]) stack.getPixels(i);
   for (int j=0;j<dimension;j++) {
    if(pixels[j] != 0.0F) {
    sum[i-1] += (float) pixels[j];
    count[i-1]++;
    }
   }
  }

  for (int i=1;i<=stack.getSize();i++) {
   average[i-1] = (float) ( sum[i-1] / count[i-1] );
   String result = "Average at " + i + " is " + average[i-1];
   result += " # of pixels: " + count[i-1];
   IJ.log(result);
  }
 }

}
