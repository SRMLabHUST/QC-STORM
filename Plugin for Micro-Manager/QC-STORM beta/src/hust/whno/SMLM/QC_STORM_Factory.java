/*
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU LESSER GENERAL PUBLIC LICENSE as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU LESSER GENERAL PUBLIC LICENSE for more details.

You should have received a copy of the GNU LESSER GENERAL PUBLIC LICENSE
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/


package hust.whno.SMLM;


import org.micromanager.data.Processor;
import org.micromanager.data.ProcessorFactory;

import org.micromanager.Studio;



public class QC_STORM_Factory implements ProcessorFactory {
   private Studio studio_;
   public QC_STORM_Configurator MyConfigurator;
   public Processor MyProcessor;

   QC_STORM_Factory (Studio studio, QC_STORM_Configurator iConfigurator)
   {
//       JOptionPane.showMessageDialog(null, "create Factory", "Hello world!",      JOptionPane.PLAIN_MESSAGE);
       studio_ = studio;
       MyConfigurator = iConfigurator;
       
   }
   
    @Override
    public Processor createProcessor() {
        // get para from PropertyMap settings and create processing based on those parameters
        MyProcessor = new QC_STORM_Processor(studio_, MyConfigurator, "");
        return MyProcessor;
    }
    
}
