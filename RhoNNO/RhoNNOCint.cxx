// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME RhoNNOCint

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "VNeuralNetPlotter.h"
#include "VNeuralNet.h"
#include "VSupervisedNet.h"
#include "VUnsupervisedNet.h"
#include "TFD.h"
#include "TSGCS.h"
#include "TSGNG.h"
#include "TPerceptron.h"
#include "TMLP.h"
#include "TXMLP.h"
#include "TNeuralNetCell.h"
#include "TGCS.h"
#include "TGNG.h"
#include "TLVQ.h"
#include "TDataServe.h"
#include "TNNK.h"
#include "TGNGTracker.h"
#include "TRadon.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void delete_VNeuralNetPlotter(void *p);
   static void deleteArray_VNeuralNetPlotter(void *p);
   static void destruct_VNeuralNetPlotter(void *p);
   static void streamer_VNeuralNetPlotter(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::VNeuralNetPlotter*)
   {
      ::VNeuralNetPlotter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::VNeuralNetPlotter >(0);
      static ::ROOT::TGenericClassInfo 
         instance("VNeuralNetPlotter", ::VNeuralNetPlotter::Class_Version(), "VNeuralNetPlotter.h", 21,
                  typeid(::VNeuralNetPlotter), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::VNeuralNetPlotter::Dictionary, isa_proxy, 16,
                  sizeof(::VNeuralNetPlotter) );
      instance.SetDelete(&delete_VNeuralNetPlotter);
      instance.SetDeleteArray(&deleteArray_VNeuralNetPlotter);
      instance.SetDestructor(&destruct_VNeuralNetPlotter);
      instance.SetStreamerFunc(&streamer_VNeuralNetPlotter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::VNeuralNetPlotter*)
   {
      return GenerateInitInstanceLocal((::VNeuralNetPlotter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::VNeuralNetPlotter*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_TSimpleNeuralNetPlotter(void *p);
   static void deleteArray_TSimpleNeuralNetPlotter(void *p);
   static void destruct_TSimpleNeuralNetPlotter(void *p);
   static void streamer_TSimpleNeuralNetPlotter(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TSimpleNeuralNetPlotter*)
   {
      ::TSimpleNeuralNetPlotter *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TSimpleNeuralNetPlotter >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TSimpleNeuralNetPlotter", ::TSimpleNeuralNetPlotter::Class_Version(), "VNeuralNetPlotter.h", 38,
                  typeid(::TSimpleNeuralNetPlotter), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TSimpleNeuralNetPlotter::Dictionary, isa_proxy, 16,
                  sizeof(::TSimpleNeuralNetPlotter) );
      instance.SetDelete(&delete_TSimpleNeuralNetPlotter);
      instance.SetDeleteArray(&deleteArray_TSimpleNeuralNetPlotter);
      instance.SetDestructor(&destruct_TSimpleNeuralNetPlotter);
      instance.SetStreamerFunc(&streamer_TSimpleNeuralNetPlotter);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TSimpleNeuralNetPlotter*)
   {
      return GenerateInitInstanceLocal((::TSimpleNeuralNetPlotter*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TSimpleNeuralNetPlotter*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TNeuralNetParameters(void *p = 0);
   static void *newArray_TNeuralNetParameters(Long_t size, void *p);
   static void delete_TNeuralNetParameters(void *p);
   static void deleteArray_TNeuralNetParameters(void *p);
   static void destruct_TNeuralNetParameters(void *p);
   static void streamer_TNeuralNetParameters(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TNeuralNetParameters*)
   {
      ::TNeuralNetParameters *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TNeuralNetParameters >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TNeuralNetParameters", ::TNeuralNetParameters::Class_Version(), "VNeuralNet.h", 27,
                  typeid(::TNeuralNetParameters), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TNeuralNetParameters::Dictionary, isa_proxy, 16,
                  sizeof(::TNeuralNetParameters) );
      instance.SetNew(&new_TNeuralNetParameters);
      instance.SetNewArray(&newArray_TNeuralNetParameters);
      instance.SetDelete(&delete_TNeuralNetParameters);
      instance.SetDeleteArray(&deleteArray_TNeuralNetParameters);
      instance.SetDestructor(&destruct_TNeuralNetParameters);
      instance.SetStreamerFunc(&streamer_TNeuralNetParameters);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TNeuralNetParameters*)
   {
      return GenerateInitInstanceLocal((::TNeuralNetParameters*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TNeuralNetParameters*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_VNeuralNet(void *p);
   static void deleteArray_VNeuralNet(void *p);
   static void destruct_VNeuralNet(void *p);
   static void streamer_VNeuralNet(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::VNeuralNet*)
   {
      ::VNeuralNet *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::VNeuralNet >(0);
      static ::ROOT::TGenericClassInfo 
         instance("VNeuralNet", ::VNeuralNet::Class_Version(), "VNeuralNet.h", 51,
                  typeid(::VNeuralNet), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::VNeuralNet::Dictionary, isa_proxy, 16,
                  sizeof(::VNeuralNet) );
      instance.SetDelete(&delete_VNeuralNet);
      instance.SetDeleteArray(&deleteArray_VNeuralNet);
      instance.SetDestructor(&destruct_VNeuralNet);
      instance.SetStreamerFunc(&streamer_VNeuralNet);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::VNeuralNet*)
   {
      return GenerateInitInstanceLocal((::VNeuralNet*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::VNeuralNet*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_VSupervisedNet(void *p);
   static void deleteArray_VSupervisedNet(void *p);
   static void destruct_VSupervisedNet(void *p);
   static void streamer_VSupervisedNet(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::VSupervisedNet*)
   {
      ::VSupervisedNet *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::VSupervisedNet >(0);
      static ::ROOT::TGenericClassInfo 
         instance("VSupervisedNet", ::VSupervisedNet::Class_Version(), "VSupervisedNet.h", 18,
                  typeid(::VSupervisedNet), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::VSupervisedNet::Dictionary, isa_proxy, 16,
                  sizeof(::VSupervisedNet) );
      instance.SetDelete(&delete_VSupervisedNet);
      instance.SetDeleteArray(&deleteArray_VSupervisedNet);
      instance.SetDestructor(&destruct_VSupervisedNet);
      instance.SetStreamerFunc(&streamer_VSupervisedNet);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::VSupervisedNet*)
   {
      return GenerateInitInstanceLocal((::VSupervisedNet*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::VSupervisedNet*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TNeuralNetCellParameters(void *p = 0);
   static void *newArray_TNeuralNetCellParameters(Long_t size, void *p);
   static void delete_TNeuralNetCellParameters(void *p);
   static void deleteArray_TNeuralNetCellParameters(void *p);
   static void destruct_TNeuralNetCellParameters(void *p);
   static void streamer_TNeuralNetCellParameters(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TNeuralNetCellParameters*)
   {
      ::TNeuralNetCellParameters *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TNeuralNetCellParameters >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TNeuralNetCellParameters", ::TNeuralNetCellParameters::Class_Version(), "TNeuralNetCell.h", 15,
                  typeid(::TNeuralNetCellParameters), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TNeuralNetCellParameters::Dictionary, isa_proxy, 16,
                  sizeof(::TNeuralNetCellParameters) );
      instance.SetNew(&new_TNeuralNetCellParameters);
      instance.SetNewArray(&newArray_TNeuralNetCellParameters);
      instance.SetDelete(&delete_TNeuralNetCellParameters);
      instance.SetDeleteArray(&deleteArray_TNeuralNetCellParameters);
      instance.SetDestructor(&destruct_TNeuralNetCellParameters);
      instance.SetStreamerFunc(&streamer_TNeuralNetCellParameters);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TNeuralNetCellParameters*)
   {
      return GenerateInitInstanceLocal((::TNeuralNetCellParameters*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TNeuralNetCellParameters*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TNeuralNetCell(void *p = 0);
   static void *newArray_TNeuralNetCell(Long_t size, void *p);
   static void delete_TNeuralNetCell(void *p);
   static void deleteArray_TNeuralNetCell(void *p);
   static void destruct_TNeuralNetCell(void *p);
   static void streamer_TNeuralNetCell(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TNeuralNetCell*)
   {
      ::TNeuralNetCell *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TNeuralNetCell >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TNeuralNetCell", ::TNeuralNetCell::Class_Version(), "TNeuralNetCell.h", 58,
                  typeid(::TNeuralNetCell), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TNeuralNetCell::Dictionary, isa_proxy, 16,
                  sizeof(::TNeuralNetCell) );
      instance.SetNew(&new_TNeuralNetCell);
      instance.SetNewArray(&newArray_TNeuralNetCell);
      instance.SetDelete(&delete_TNeuralNetCell);
      instance.SetDeleteArray(&deleteArray_TNeuralNetCell);
      instance.SetDestructor(&destruct_TNeuralNetCell);
      instance.SetStreamerFunc(&streamer_TNeuralNetCell);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TNeuralNetCell*)
   {
      return GenerateInitInstanceLocal((::TNeuralNetCell*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TNeuralNetCell*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void delete_VUnsupervisedNet(void *p);
   static void deleteArray_VUnsupervisedNet(void *p);
   static void destruct_VUnsupervisedNet(void *p);
   static void streamer_VUnsupervisedNet(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::VUnsupervisedNet*)
   {
      ::VUnsupervisedNet *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::VUnsupervisedNet >(0);
      static ::ROOT::TGenericClassInfo 
         instance("VUnsupervisedNet", ::VUnsupervisedNet::Class_Version(), "VUnsupervisedNet.h", 19,
                  typeid(::VUnsupervisedNet), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::VUnsupervisedNet::Dictionary, isa_proxy, 16,
                  sizeof(::VUnsupervisedNet) );
      instance.SetDelete(&delete_VUnsupervisedNet);
      instance.SetDeleteArray(&deleteArray_VUnsupervisedNet);
      instance.SetDestructor(&destruct_VUnsupervisedNet);
      instance.SetStreamerFunc(&streamer_VUnsupervisedNet);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::VUnsupervisedNet*)
   {
      return GenerateInitInstanceLocal((::VUnsupervisedNet*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::VUnsupervisedNet*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TPerceptron(void *p = 0);
   static void *newArray_TPerceptron(Long_t size, void *p);
   static void delete_TPerceptron(void *p);
   static void deleteArray_TPerceptron(void *p);
   static void destruct_TPerceptron(void *p);
   static void streamer_TPerceptron(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TPerceptron*)
   {
      ::TPerceptron *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TPerceptron >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TPerceptron", ::TPerceptron::Class_Version(), "TPerceptron.h", 27,
                  typeid(::TPerceptron), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TPerceptron::Dictionary, isa_proxy, 16,
                  sizeof(::TPerceptron) );
      instance.SetNew(&new_TPerceptron);
      instance.SetNewArray(&newArray_TPerceptron);
      instance.SetDelete(&delete_TPerceptron);
      instance.SetDeleteArray(&deleteArray_TPerceptron);
      instance.SetDestructor(&destruct_TPerceptron);
      instance.SetStreamerFunc(&streamer_TPerceptron);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TPerceptron*)
   {
      return GenerateInitInstanceLocal((::TPerceptron*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TPerceptron*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TXMLP(void *p = 0);
   static void *newArray_TXMLP(Long_t size, void *p);
   static void delete_TXMLP(void *p);
   static void deleteArray_TXMLP(void *p);
   static void destruct_TXMLP(void *p);
   static void streamer_TXMLP(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TXMLP*)
   {
      ::TXMLP *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TXMLP >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TXMLP", ::TXMLP::Class_Version(), "TXMLP.h", 16,
                  typeid(::TXMLP), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TXMLP::Dictionary, isa_proxy, 16,
                  sizeof(::TXMLP) );
      instance.SetNew(&new_TXMLP);
      instance.SetNewArray(&newArray_TXMLP);
      instance.SetDelete(&delete_TXMLP);
      instance.SetDeleteArray(&deleteArray_TXMLP);
      instance.SetDestructor(&destruct_TXMLP);
      instance.SetStreamerFunc(&streamer_TXMLP);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TXMLP*)
   {
      return GenerateInitInstanceLocal((::TXMLP*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TXMLP*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TFD(void *p = 0);
   static void *newArray_TFD(Long_t size, void *p);
   static void delete_TFD(void *p);
   static void deleteArray_TFD(void *p);
   static void destruct_TFD(void *p);
   static void streamer_TFD(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TFD*)
   {
      ::TFD *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TFD >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TFD", ::TFD::Class_Version(), "TFD.h", 15,
                  typeid(::TFD), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TFD::Dictionary, isa_proxy, 16,
                  sizeof(::TFD) );
      instance.SetNew(&new_TFD);
      instance.SetNewArray(&newArray_TFD);
      instance.SetDelete(&delete_TFD);
      instance.SetDeleteArray(&deleteArray_TFD);
      instance.SetDestructor(&destruct_TFD);
      instance.SetStreamerFunc(&streamer_TFD);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TFD*)
   {
      return GenerateInitInstanceLocal((::TFD*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TFD*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TSGCS(void *p = 0);
   static void *newArray_TSGCS(Long_t size, void *p);
   static void delete_TSGCS(void *p);
   static void deleteArray_TSGCS(void *p);
   static void destruct_TSGCS(void *p);
   static void streamer_TSGCS(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TSGCS*)
   {
      ::TSGCS *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TSGCS >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TSGCS", ::TSGCS::Class_Version(), "TSGCS.h", 16,
                  typeid(::TSGCS), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TSGCS::Dictionary, isa_proxy, 16,
                  sizeof(::TSGCS) );
      instance.SetNew(&new_TSGCS);
      instance.SetNewArray(&newArray_TSGCS);
      instance.SetDelete(&delete_TSGCS);
      instance.SetDeleteArray(&deleteArray_TSGCS);
      instance.SetDestructor(&destruct_TSGCS);
      instance.SetStreamerFunc(&streamer_TSGCS);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TSGCS*)
   {
      return GenerateInitInstanceLocal((::TSGCS*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TSGCS*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TSGNG(void *p = 0);
   static void *newArray_TSGNG(Long_t size, void *p);
   static void delete_TSGNG(void *p);
   static void deleteArray_TSGNG(void *p);
   static void destruct_TSGNG(void *p);
   static void streamer_TSGNG(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TSGNG*)
   {
      ::TSGNG *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TSGNG >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TSGNG", ::TSGNG::Class_Version(), "TSGNG.h", 16,
                  typeid(::TSGNG), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TSGNG::Dictionary, isa_proxy, 16,
                  sizeof(::TSGNG) );
      instance.SetNew(&new_TSGNG);
      instance.SetNewArray(&newArray_TSGNG);
      instance.SetDelete(&delete_TSGNG);
      instance.SetDeleteArray(&deleteArray_TSGNG);
      instance.SetDestructor(&destruct_TSGNG);
      instance.SetStreamerFunc(&streamer_TSGNG);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TSGNG*)
   {
      return GenerateInitInstanceLocal((::TSGNG*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TSGNG*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TMLP(void *p = 0);
   static void *newArray_TMLP(Long_t size, void *p);
   static void delete_TMLP(void *p);
   static void deleteArray_TMLP(void *p);
   static void destruct_TMLP(void *p);
   static void streamer_TMLP(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TMLP*)
   {
      ::TMLP *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TMLP >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TMLP", ::TMLP::Class_Version(), "TMLP.h", 15,
                  typeid(::TMLP), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TMLP::Dictionary, isa_proxy, 16,
                  sizeof(::TMLP) );
      instance.SetNew(&new_TMLP);
      instance.SetNewArray(&newArray_TMLP);
      instance.SetDelete(&delete_TMLP);
      instance.SetDeleteArray(&deleteArray_TMLP);
      instance.SetDestructor(&destruct_TMLP);
      instance.SetStreamerFunc(&streamer_TMLP);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TMLP*)
   {
      return GenerateInitInstanceLocal((::TMLP*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TMLP*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TGCS(void *p = 0);
   static void *newArray_TGCS(Long_t size, void *p);
   static void delete_TGCS(void *p);
   static void deleteArray_TGCS(void *p);
   static void destruct_TGCS(void *p);
   static void streamer_TGCS(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TGCS*)
   {
      ::TGCS *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TGCS >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TGCS", ::TGCS::Class_Version(), "TGCS.h", 15,
                  typeid(::TGCS), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TGCS::Dictionary, isa_proxy, 16,
                  sizeof(::TGCS) );
      instance.SetNew(&new_TGCS);
      instance.SetNewArray(&newArray_TGCS);
      instance.SetDelete(&delete_TGCS);
      instance.SetDeleteArray(&deleteArray_TGCS);
      instance.SetDestructor(&destruct_TGCS);
      instance.SetStreamerFunc(&streamer_TGCS);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TGCS*)
   {
      return GenerateInitInstanceLocal((::TGCS*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TGCS*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TGNG(void *p = 0);
   static void *newArray_TGNG(Long_t size, void *p);
   static void delete_TGNG(void *p);
   static void deleteArray_TGNG(void *p);
   static void destruct_TGNG(void *p);
   static void streamer_TGNG(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TGNG*)
   {
      ::TGNG *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TGNG >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TGNG", ::TGNG::Class_Version(), "TGNG.h", 15,
                  typeid(::TGNG), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TGNG::Dictionary, isa_proxy, 16,
                  sizeof(::TGNG) );
      instance.SetNew(&new_TGNG);
      instance.SetNewArray(&newArray_TGNG);
      instance.SetDelete(&delete_TGNG);
      instance.SetDeleteArray(&deleteArray_TGNG);
      instance.SetDestructor(&destruct_TGNG);
      instance.SetStreamerFunc(&streamer_TGNG);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TGNG*)
   {
      return GenerateInitInstanceLocal((::TGNG*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TGNG*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TLVQ(void *p = 0);
   static void *newArray_TLVQ(Long_t size, void *p);
   static void delete_TLVQ(void *p);
   static void deleteArray_TLVQ(void *p);
   static void destruct_TLVQ(void *p);
   static void streamer_TLVQ(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TLVQ*)
   {
      ::TLVQ *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TLVQ >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TLVQ", ::TLVQ::Class_Version(), "TLVQ.h", 15,
                  typeid(::TLVQ), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TLVQ::Dictionary, isa_proxy, 16,
                  sizeof(::TLVQ) );
      instance.SetNew(&new_TLVQ);
      instance.SetNewArray(&newArray_TLVQ);
      instance.SetDelete(&delete_TLVQ);
      instance.SetDeleteArray(&deleteArray_TLVQ);
      instance.SetDestructor(&destruct_TLVQ);
      instance.SetStreamerFunc(&streamer_TLVQ);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TLVQ*)
   {
      return GenerateInitInstanceLocal((::TLVQ*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TLVQ*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TDataServe(void *p = 0);
   static void *newArray_TDataServe(Long_t size, void *p);
   static void delete_TDataServe(void *p);
   static void deleteArray_TDataServe(void *p);
   static void destruct_TDataServe(void *p);
   static void streamer_TDataServe(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TDataServe*)
   {
      ::TDataServe *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TDataServe >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TDataServe", ::TDataServe::Class_Version(), "TDataServe.h", 30,
                  typeid(::TDataServe), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TDataServe::Dictionary, isa_proxy, 17,
                  sizeof(::TDataServe) );
      instance.SetNew(&new_TDataServe);
      instance.SetNewArray(&newArray_TDataServe);
      instance.SetDelete(&delete_TDataServe);
      instance.SetDeleteArray(&deleteArray_TDataServe);
      instance.SetDestructor(&destruct_TDataServe);
      instance.SetStreamerFunc(&streamer_TDataServe);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TDataServe*)
   {
      return GenerateInitInstanceLocal((::TDataServe*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TDataServe*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TNNK(void *p = 0);
   static void *newArray_TNNK(Long_t size, void *p);
   static void delete_TNNK(void *p);
   static void deleteArray_TNNK(void *p);
   static void destruct_TNNK(void *p);
   static void streamer_TNNK(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TNNK*)
   {
      ::TNNK *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TNNK >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TNNK", ::TNNK::Class_Version(), "TNNK.h", 17,
                  typeid(::TNNK), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TNNK::Dictionary, isa_proxy, 16,
                  sizeof(::TNNK) );
      instance.SetNew(&new_TNNK);
      instance.SetNewArray(&newArray_TNNK);
      instance.SetDelete(&delete_TNNK);
      instance.SetDeleteArray(&deleteArray_TNNK);
      instance.SetDestructor(&destruct_TNNK);
      instance.SetStreamerFunc(&streamer_TNNK);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TNNK*)
   {
      return GenerateInitInstanceLocal((::TNNK*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TNNK*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TNNKernel(void *p = 0);
   static void *newArray_TNNKernel(Long_t size, void *p);
   static void delete_TNNKernel(void *p);
   static void deleteArray_TNNKernel(void *p);
   static void destruct_TNNKernel(void *p);
   static void streamer_TNNKernel(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TNNKernel*)
   {
      ::TNNKernel *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TNNKernel >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TNNKernel", ::TNNKernel::Class_Version(), "TNNK.h", 53,
                  typeid(::TNNKernel), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TNNKernel::Dictionary, isa_proxy, 16,
                  sizeof(::TNNKernel) );
      instance.SetNew(&new_TNNKernel);
      instance.SetNewArray(&newArray_TNNKernel);
      instance.SetDelete(&delete_TNNKernel);
      instance.SetDeleteArray(&deleteArray_TNNKernel);
      instance.SetDestructor(&destruct_TNNKernel);
      instance.SetStreamerFunc(&streamer_TNNKernel);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TNNKernel*)
   {
      return GenerateInitInstanceLocal((::TNNKernel*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TNNKernel*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TGNGTracker(void *p = 0);
   static void *newArray_TGNGTracker(Long_t size, void *p);
   static void delete_TGNGTracker(void *p);
   static void deleteArray_TGNGTracker(void *p);
   static void destruct_TGNGTracker(void *p);
   static void streamer_TGNGTracker(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TGNGTracker*)
   {
      ::TGNGTracker *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TGNGTracker >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TGNGTracker", ::TGNGTracker::Class_Version(), "TGNGTracker.h", 17,
                  typeid(::TGNGTracker), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TGNGTracker::Dictionary, isa_proxy, 16,
                  sizeof(::TGNGTracker) );
      instance.SetNew(&new_TGNGTracker);
      instance.SetNewArray(&newArray_TGNGTracker);
      instance.SetDelete(&delete_TGNGTracker);
      instance.SetDeleteArray(&deleteArray_TGNGTracker);
      instance.SetDestructor(&destruct_TGNGTracker);
      instance.SetStreamerFunc(&streamer_TGNGTracker);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TGNGTracker*)
   {
      return GenerateInitInstanceLocal((::TGNGTracker*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TGNGTracker*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TRadon(void *p = 0);
   static void *newArray_TRadon(Long_t size, void *p);
   static void delete_TRadon(void *p);
   static void deleteArray_TRadon(void *p);
   static void destruct_TRadon(void *p);
   static void streamer_TRadon(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TRadon*)
   {
      ::TRadon *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TRadon >(0);
      static ::ROOT::TGenericClassInfo 
         instance("TRadon", ::TRadon::Class_Version(), "TRadon.h", 21,
                  typeid(::TRadon), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TRadon::Dictionary, isa_proxy, 16,
                  sizeof(::TRadon) );
      instance.SetNew(&new_TRadon);
      instance.SetNewArray(&newArray_TRadon);
      instance.SetDelete(&delete_TRadon);
      instance.SetDeleteArray(&deleteArray_TRadon);
      instance.SetDestructor(&destruct_TRadon);
      instance.SetStreamerFunc(&streamer_TRadon);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TRadon*)
   {
      return GenerateInitInstanceLocal((::TRadon*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::TRadon*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr VNeuralNetPlotter::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *VNeuralNetPlotter::Class_Name()
{
   return "VNeuralNetPlotter";
}

//______________________________________________________________________________
const char *VNeuralNetPlotter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::VNeuralNetPlotter*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int VNeuralNetPlotter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::VNeuralNetPlotter*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *VNeuralNetPlotter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::VNeuralNetPlotter*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *VNeuralNetPlotter::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::VNeuralNetPlotter*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TSimpleNeuralNetPlotter::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TSimpleNeuralNetPlotter::Class_Name()
{
   return "TSimpleNeuralNetPlotter";
}

//______________________________________________________________________________
const char *TSimpleNeuralNetPlotter::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TSimpleNeuralNetPlotter*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TSimpleNeuralNetPlotter::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TSimpleNeuralNetPlotter*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TSimpleNeuralNetPlotter::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TSimpleNeuralNetPlotter*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TSimpleNeuralNetPlotter::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TSimpleNeuralNetPlotter*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TNeuralNetParameters::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TNeuralNetParameters::Class_Name()
{
   return "TNeuralNetParameters";
}

//______________________________________________________________________________
const char *TNeuralNetParameters::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TNeuralNetParameters*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TNeuralNetParameters::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TNeuralNetParameters*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TNeuralNetParameters::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TNeuralNetParameters*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TNeuralNetParameters::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TNeuralNetParameters*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr VNeuralNet::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *VNeuralNet::Class_Name()
{
   return "VNeuralNet";
}

//______________________________________________________________________________
const char *VNeuralNet::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::VNeuralNet*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int VNeuralNet::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::VNeuralNet*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *VNeuralNet::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::VNeuralNet*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *VNeuralNet::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::VNeuralNet*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr VSupervisedNet::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *VSupervisedNet::Class_Name()
{
   return "VSupervisedNet";
}

//______________________________________________________________________________
const char *VSupervisedNet::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::VSupervisedNet*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int VSupervisedNet::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::VSupervisedNet*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *VSupervisedNet::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::VSupervisedNet*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *VSupervisedNet::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::VSupervisedNet*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TNeuralNetCellParameters::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TNeuralNetCellParameters::Class_Name()
{
   return "TNeuralNetCellParameters";
}

//______________________________________________________________________________
const char *TNeuralNetCellParameters::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TNeuralNetCellParameters*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TNeuralNetCellParameters::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TNeuralNetCellParameters*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TNeuralNetCellParameters::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TNeuralNetCellParameters*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TNeuralNetCellParameters::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TNeuralNetCellParameters*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TNeuralNetCell::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TNeuralNetCell::Class_Name()
{
   return "TNeuralNetCell";
}

//______________________________________________________________________________
const char *TNeuralNetCell::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TNeuralNetCell*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TNeuralNetCell::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TNeuralNetCell*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TNeuralNetCell::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TNeuralNetCell*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TNeuralNetCell::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TNeuralNetCell*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr VUnsupervisedNet::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *VUnsupervisedNet::Class_Name()
{
   return "VUnsupervisedNet";
}

//______________________________________________________________________________
const char *VUnsupervisedNet::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::VUnsupervisedNet*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int VUnsupervisedNet::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::VUnsupervisedNet*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *VUnsupervisedNet::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::VUnsupervisedNet*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *VUnsupervisedNet::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::VUnsupervisedNet*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TPerceptron::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TPerceptron::Class_Name()
{
   return "TPerceptron";
}

//______________________________________________________________________________
const char *TPerceptron::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TPerceptron*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TPerceptron::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TPerceptron*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TPerceptron::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TPerceptron*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TPerceptron::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TPerceptron*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TXMLP::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TXMLP::Class_Name()
{
   return "TXMLP";
}

//______________________________________________________________________________
const char *TXMLP::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TXMLP*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TXMLP::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TXMLP*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TXMLP::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TXMLP*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TXMLP::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TXMLP*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TFD::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TFD::Class_Name()
{
   return "TFD";
}

//______________________________________________________________________________
const char *TFD::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TFD*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TFD::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TFD*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TFD::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TFD*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TFD::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TFD*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TSGCS::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TSGCS::Class_Name()
{
   return "TSGCS";
}

//______________________________________________________________________________
const char *TSGCS::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TSGCS*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TSGCS::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TSGCS*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TSGCS::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TSGCS*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TSGCS::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TSGCS*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TSGNG::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TSGNG::Class_Name()
{
   return "TSGNG";
}

//______________________________________________________________________________
const char *TSGNG::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TSGNG*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TSGNG::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TSGNG*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TSGNG::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TSGNG*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TSGNG::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TSGNG*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TMLP::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TMLP::Class_Name()
{
   return "TMLP";
}

//______________________________________________________________________________
const char *TMLP::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TMLP*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TMLP::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TMLP*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TMLP::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TMLP*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TMLP::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TMLP*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TGCS::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TGCS::Class_Name()
{
   return "TGCS";
}

//______________________________________________________________________________
const char *TGCS::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TGCS*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TGCS::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TGCS*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TGCS::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TGCS*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TGCS::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TGCS*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TGNG::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TGNG::Class_Name()
{
   return "TGNG";
}

//______________________________________________________________________________
const char *TGNG::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TGNG*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TGNG::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TGNG*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TGNG::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TGNG*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TGNG::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TGNG*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TLVQ::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TLVQ::Class_Name()
{
   return "TLVQ";
}

//______________________________________________________________________________
const char *TLVQ::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TLVQ*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TLVQ::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TLVQ*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TLVQ::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TLVQ*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TLVQ::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TLVQ*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TDataServe::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TDataServe::Class_Name()
{
   return "TDataServe";
}

//______________________________________________________________________________
const char *TDataServe::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TDataServe*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TDataServe::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TDataServe*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TDataServe::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TDataServe*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TDataServe::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TDataServe*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TNNK::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TNNK::Class_Name()
{
   return "TNNK";
}

//______________________________________________________________________________
const char *TNNK::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TNNK*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TNNK::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TNNK*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TNNK::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TNNK*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TNNK::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TNNK*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TNNKernel::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TNNKernel::Class_Name()
{
   return "TNNKernel";
}

//______________________________________________________________________________
const char *TNNKernel::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TNNKernel*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TNNKernel::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TNNKernel*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TNNKernel::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TNNKernel*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TNNKernel::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TNNKernel*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TGNGTracker::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TGNGTracker::Class_Name()
{
   return "TGNGTracker";
}

//______________________________________________________________________________
const char *TGNGTracker::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TGNGTracker*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TGNGTracker::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TGNGTracker*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TGNGTracker::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TGNGTracker*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TGNGTracker::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TGNGTracker*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TRadon::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *TRadon::Class_Name()
{
   return "TRadon";
}

//______________________________________________________________________________
const char *TRadon::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRadon*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int TRadon::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TRadon*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TRadon::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRadon*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TRadon::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TRadon*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void VNeuralNetPlotter::Streamer(TBuffer &R__b)
{
   // Stream an object of class VNeuralNetPlotter.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, VNeuralNetPlotter::IsA());
   } else {
      R__c = R__b.WriteVersion(VNeuralNetPlotter::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_VNeuralNetPlotter(void *p) {
      delete ((::VNeuralNetPlotter*)p);
   }
   static void deleteArray_VNeuralNetPlotter(void *p) {
      delete [] ((::VNeuralNetPlotter*)p);
   }
   static void destruct_VNeuralNetPlotter(void *p) {
      typedef ::VNeuralNetPlotter current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_VNeuralNetPlotter(TBuffer &buf, void *obj) {
      ((::VNeuralNetPlotter*)obj)->::VNeuralNetPlotter::Streamer(buf);
   }
} // end of namespace ROOT for class ::VNeuralNetPlotter

//______________________________________________________________________________
void TSimpleNeuralNetPlotter::Streamer(TBuffer &R__b)
{
   // Stream an object of class TSimpleNeuralNetPlotter.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      VNeuralNetPlotter::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, TSimpleNeuralNetPlotter::IsA());
   } else {
      R__c = R__b.WriteVersion(TSimpleNeuralNetPlotter::IsA(), kTRUE);
      VNeuralNetPlotter::Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_TSimpleNeuralNetPlotter(void *p) {
      delete ((::TSimpleNeuralNetPlotter*)p);
   }
   static void deleteArray_TSimpleNeuralNetPlotter(void *p) {
      delete [] ((::TSimpleNeuralNetPlotter*)p);
   }
   static void destruct_TSimpleNeuralNetPlotter(void *p) {
      typedef ::TSimpleNeuralNetPlotter current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TSimpleNeuralNetPlotter(TBuffer &buf, void *obj) {
      ((::TSimpleNeuralNetPlotter*)obj)->::TSimpleNeuralNetPlotter::Streamer(buf);
   }
} // end of namespace ROOT for class ::TSimpleNeuralNetPlotter

//______________________________________________________________________________
void TNeuralNetParameters::Streamer(TBuffer &R__b)
{
   // Stream an object of class TNeuralNetParameters.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b.ReadStaticArray((char*)fNetId);
      R__b >> fLayers;
      R__b >> fInScale;
      R__b >> fInNodes;
      R__b >> fOutNodes;
      R__b >> fLearnStep;
      R__b >> fMu;
      R__b >> fFse;
      void *ptr_fTransferId = (void*)&fTransferId;
      R__b >> *reinterpret_cast<Int_t*>(ptr_fTransferId);
      R__b >> fPerceptronId;
      R__b >> fThreshold;
      R__b.CheckByteCount(R__s, R__c, TNeuralNetParameters::IsA());
   } else {
      R__c = R__b.WriteVersion(TNeuralNetParameters::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b.WriteArray(fNetId, 9);
      R__b << fLayers;
      R__b << fInScale;
      R__b << fInNodes;
      R__b << fOutNodes;
      R__b << fLearnStep;
      R__b << fMu;
      R__b << fFse;
      R__b << (Int_t)fTransferId;
      R__b << fPerceptronId;
      R__b << fThreshold;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TNeuralNetParameters(void *p) {
      return  p ? new(p) ::TNeuralNetParameters : new ::TNeuralNetParameters;
   }
   static void *newArray_TNeuralNetParameters(Long_t nElements, void *p) {
      return p ? new(p) ::TNeuralNetParameters[nElements] : new ::TNeuralNetParameters[nElements];
   }
   // Wrapper around operator delete
   static void delete_TNeuralNetParameters(void *p) {
      delete ((::TNeuralNetParameters*)p);
   }
   static void deleteArray_TNeuralNetParameters(void *p) {
      delete [] ((::TNeuralNetParameters*)p);
   }
   static void destruct_TNeuralNetParameters(void *p) {
      typedef ::TNeuralNetParameters current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TNeuralNetParameters(TBuffer &buf, void *obj) {
      ((::TNeuralNetParameters*)obj)->::TNeuralNetParameters::Streamer(buf);
   }
} // end of namespace ROOT for class ::TNeuralNetParameters

//______________________________________________________________________________
void VNeuralNet::Streamer(TBuffer &R__b)
{
   // Stream an object of class VNeuralNet.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      { TString R__str; R__str.Streamer(R__b); fFilename = R__str.Data(); }
      void *ptr_fFiletype = (void*)&fFiletype;
      R__b >> *reinterpret_cast<Int_t*>(ptr_fFiletype);
      fParm.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, VNeuralNet::IsA());
   } else {
      R__c = R__b.WriteVersion(VNeuralNet::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      { TString R__str = fFilename.c_str(); R__str.Streamer(R__b);}
      R__b << (Int_t)fFiletype;
      fParm.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_VNeuralNet(void *p) {
      delete ((::VNeuralNet*)p);
   }
   static void deleteArray_VNeuralNet(void *p) {
      delete [] ((::VNeuralNet*)p);
   }
   static void destruct_VNeuralNet(void *p) {
      typedef ::VNeuralNet current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_VNeuralNet(TBuffer &buf, void *obj) {
      ((::VNeuralNet*)obj)->::VNeuralNet::Streamer(buf);
   }
} // end of namespace ROOT for class ::VNeuralNet

//______________________________________________________________________________
void VSupervisedNet::Streamer(TBuffer &R__b)
{
   // Stream an object of class VSupervisedNet.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      VNeuralNet::Streamer(R__b);
      R__b >> fTuple;
      R__b.CheckByteCount(R__s, R__c, VSupervisedNet::IsA());
   } else {
      R__c = R__b.WriteVersion(VSupervisedNet::IsA(), kTRUE);
      VNeuralNet::Streamer(R__b);
      R__b << fTuple;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_VSupervisedNet(void *p) {
      delete ((::VSupervisedNet*)p);
   }
   static void deleteArray_VSupervisedNet(void *p) {
      delete [] ((::VSupervisedNet*)p);
   }
   static void destruct_VSupervisedNet(void *p) {
      typedef ::VSupervisedNet current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_VSupervisedNet(TBuffer &buf, void *obj) {
      ((::VSupervisedNet*)obj)->::VSupervisedNet::Streamer(buf);
   }
} // end of namespace ROOT for class ::VSupervisedNet

//______________________________________________________________________________
void TNeuralNetCellParameters::Streamer(TBuffer &R__b)
{
   // Stream an object of class TNeuralNetCellParameters.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> fWinStep;
      R__b >> fNeiStep;
      R__b >> fNeuStep;
      R__b >> fCells;
      R__b >> fConnectors;
      R__b >> fWinCount;
      R__b >> fMinCells;
      R__b >> fMaxCells;
      R__b >> fInsertStep;
      R__b >> fDeleteStep;
      R__b >> fInsertCount;
      R__b >> fDeleteCount;
      R__b >> fEdgeCount;
      R__b >> fErrCount;
      R__b >> fNeiCount;
      R__b >> fMinCount;
      R__b >> fMainWinCount;
      R__b >> fMainErrCount;
      R__b >> fMainEdgeCount;
      R__b.CheckByteCount(R__s, R__c, TNeuralNetCellParameters::IsA());
   } else {
      R__c = R__b.WriteVersion(TNeuralNetCellParameters::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << fWinStep;
      R__b << fNeiStep;
      R__b << fNeuStep;
      R__b << fCells;
      R__b << fConnectors;
      R__b << fWinCount;
      R__b << fMinCells;
      R__b << fMaxCells;
      R__b << fInsertStep;
      R__b << fDeleteStep;
      R__b << fInsertCount;
      R__b << fDeleteCount;
      R__b << fEdgeCount;
      R__b << fErrCount;
      R__b << fNeiCount;
      R__b << fMinCount;
      R__b << fMainWinCount;
      R__b << fMainErrCount;
      R__b << fMainEdgeCount;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TNeuralNetCellParameters(void *p) {
      return  p ? new(p) ::TNeuralNetCellParameters : new ::TNeuralNetCellParameters;
   }
   static void *newArray_TNeuralNetCellParameters(Long_t nElements, void *p) {
      return p ? new(p) ::TNeuralNetCellParameters[nElements] : new ::TNeuralNetCellParameters[nElements];
   }
   // Wrapper around operator delete
   static void delete_TNeuralNetCellParameters(void *p) {
      delete ((::TNeuralNetCellParameters*)p);
   }
   static void deleteArray_TNeuralNetCellParameters(void *p) {
      delete [] ((::TNeuralNetCellParameters*)p);
   }
   static void destruct_TNeuralNetCellParameters(void *p) {
      typedef ::TNeuralNetCellParameters current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TNeuralNetCellParameters(TBuffer &buf, void *obj) {
      ((::TNeuralNetCellParameters*)obj)->::TNeuralNetCellParameters::Streamer(buf);
   }
} // end of namespace ROOT for class ::TNeuralNetCellParameters

//______________________________________________________________________________
void TNeuralNetCell::Streamer(TBuffer &R__b)
{
   // Stream an object of class TNeuralNetCell.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> fID;
      R__b >> fCount;
      R__b >> fState;
      R__b >> fNc;
      R__b >> fChi2;
      R__b >> fOut;
      R__b >> fClass;
      R__b.CheckByteCount(R__s, R__c, TNeuralNetCell::IsA());
   } else {
      R__c = R__b.WriteVersion(TNeuralNetCell::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << fID;
      R__b << fCount;
      R__b << fState;
      R__b << fNc;
      R__b << fChi2;
      R__b << fOut;
      R__b << fClass;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TNeuralNetCell(void *p) {
      return  p ? new(p) ::TNeuralNetCell : new ::TNeuralNetCell;
   }
   static void *newArray_TNeuralNetCell(Long_t nElements, void *p) {
      return p ? new(p) ::TNeuralNetCell[nElements] : new ::TNeuralNetCell[nElements];
   }
   // Wrapper around operator delete
   static void delete_TNeuralNetCell(void *p) {
      delete ((::TNeuralNetCell*)p);
   }
   static void deleteArray_TNeuralNetCell(void *p) {
      delete [] ((::TNeuralNetCell*)p);
   }
   static void destruct_TNeuralNetCell(void *p) {
      typedef ::TNeuralNetCell current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TNeuralNetCell(TBuffer &buf, void *obj) {
      ((::TNeuralNetCell*)obj)->::TNeuralNetCell::Streamer(buf);
   }
} // end of namespace ROOT for class ::TNeuralNetCell

//______________________________________________________________________________
void VUnsupervisedNet::Streamer(TBuffer &R__b)
{
   // Stream an object of class VUnsupervisedNet.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      VNeuralNet::Streamer(R__b);
      fXB.Streamer(R__b);
      R__b >> fTuple;
      R__b.CheckByteCount(R__s, R__c, VUnsupervisedNet::IsA());
   } else {
      R__c = R__b.WriteVersion(VUnsupervisedNet::IsA(), kTRUE);
      VNeuralNet::Streamer(R__b);
      fXB.Streamer(R__b);
      R__b << fTuple;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrapper around operator delete
   static void delete_VUnsupervisedNet(void *p) {
      delete ((::VUnsupervisedNet*)p);
   }
   static void deleteArray_VUnsupervisedNet(void *p) {
      delete [] ((::VUnsupervisedNet*)p);
   }
   static void destruct_VUnsupervisedNet(void *p) {
      typedef ::VUnsupervisedNet current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_VUnsupervisedNet(TBuffer &buf, void *obj) {
      ((::VUnsupervisedNet*)obj)->::VUnsupervisedNet::Streamer(buf);
   }
} // end of namespace ROOT for class ::VUnsupervisedNet

//______________________________________________________________________________
void TPerceptron::Streamer(TBuffer &R__b)
{
   // Stream an object of class TPerceptron.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      VSupervisedNet::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, TPerceptron::IsA());
   } else {
      R__c = R__b.WriteVersion(TPerceptron::IsA(), kTRUE);
      VSupervisedNet::Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TPerceptron(void *p) {
      return  p ? new(p) ::TPerceptron : new ::TPerceptron;
   }
   static void *newArray_TPerceptron(Long_t nElements, void *p) {
      return p ? new(p) ::TPerceptron[nElements] : new ::TPerceptron[nElements];
   }
   // Wrapper around operator delete
   static void delete_TPerceptron(void *p) {
      delete ((::TPerceptron*)p);
   }
   static void deleteArray_TPerceptron(void *p) {
      delete [] ((::TPerceptron*)p);
   }
   static void destruct_TPerceptron(void *p) {
      typedef ::TPerceptron current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TPerceptron(TBuffer &buf, void *obj) {
      ((::TPerceptron*)obj)->::TPerceptron::Streamer(buf);
   }
} // end of namespace ROOT for class ::TPerceptron

//______________________________________________________________________________
void TXMLP::Streamer(TBuffer &R__b)
{
   // Stream an object of class TXMLP.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      VSupervisedNet::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, TXMLP::IsA());
   } else {
      R__c = R__b.WriteVersion(TXMLP::IsA(), kTRUE);
      VSupervisedNet::Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TXMLP(void *p) {
      return  p ? new(p) ::TXMLP : new ::TXMLP;
   }
   static void *newArray_TXMLP(Long_t nElements, void *p) {
      return p ? new(p) ::TXMLP[nElements] : new ::TXMLP[nElements];
   }
   // Wrapper around operator delete
   static void delete_TXMLP(void *p) {
      delete ((::TXMLP*)p);
   }
   static void deleteArray_TXMLP(void *p) {
      delete [] ((::TXMLP*)p);
   }
   static void destruct_TXMLP(void *p) {
      typedef ::TXMLP current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TXMLP(TBuffer &buf, void *obj) {
      ((::TXMLP*)obj)->::TXMLP::Streamer(buf);
   }
} // end of namespace ROOT for class ::TXMLP

//______________________________________________________________________________
void TFD::Streamer(TBuffer &R__b)
{
   // Stream an object of class TFD.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TXMLP::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, TFD::IsA());
   } else {
      R__c = R__b.WriteVersion(TFD::IsA(), kTRUE);
      TXMLP::Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TFD(void *p) {
      return  p ? new(p) ::TFD : new ::TFD;
   }
   static void *newArray_TFD(Long_t nElements, void *p) {
      return p ? new(p) ::TFD[nElements] : new ::TFD[nElements];
   }
   // Wrapper around operator delete
   static void delete_TFD(void *p) {
      delete ((::TFD*)p);
   }
   static void deleteArray_TFD(void *p) {
      delete [] ((::TFD*)p);
   }
   static void destruct_TFD(void *p) {
      typedef ::TFD current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TFD(TBuffer &buf, void *obj) {
      ((::TFD*)obj)->::TFD::Streamer(buf);
   }
} // end of namespace ROOT for class ::TFD

//______________________________________________________________________________
void TSGCS::Streamer(TBuffer &R__b)
{
   // Stream an object of class TSGCS.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      VSupervisedNet::Streamer(R__b);
      fXB.Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, TSGCS::IsA());
   } else {
      R__c = R__b.WriteVersion(TSGCS::IsA(), kTRUE);
      VSupervisedNet::Streamer(R__b);
      fXB.Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TSGCS(void *p) {
      return  p ? new(p) ::TSGCS : new ::TSGCS;
   }
   static void *newArray_TSGCS(Long_t nElements, void *p) {
      return p ? new(p) ::TSGCS[nElements] : new ::TSGCS[nElements];
   }
   // Wrapper around operator delete
   static void delete_TSGCS(void *p) {
      delete ((::TSGCS*)p);
   }
   static void deleteArray_TSGCS(void *p) {
      delete [] ((::TSGCS*)p);
   }
   static void destruct_TSGCS(void *p) {
      typedef ::TSGCS current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TSGCS(TBuffer &buf, void *obj) {
      ((::TSGCS*)obj)->::TSGCS::Streamer(buf);
   }
} // end of namespace ROOT for class ::TSGCS

//______________________________________________________________________________
void TSGNG::Streamer(TBuffer &R__b)
{
   // Stream an object of class TSGNG.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      VSupervisedNet::Streamer(R__b);
      fXB.Streamer(R__b);
      R__b >> fMinDistSquare1;
      R__b >> fMinDistSquare2;
      R__b.CheckByteCount(R__s, R__c, TSGNG::IsA());
   } else {
      R__c = R__b.WriteVersion(TSGNG::IsA(), kTRUE);
      VSupervisedNet::Streamer(R__b);
      fXB.Streamer(R__b);
      R__b << fMinDistSquare1;
      R__b << fMinDistSquare2;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TSGNG(void *p) {
      return  p ? new(p) ::TSGNG : new ::TSGNG;
   }
   static void *newArray_TSGNG(Long_t nElements, void *p) {
      return p ? new(p) ::TSGNG[nElements] : new ::TSGNG[nElements];
   }
   // Wrapper around operator delete
   static void delete_TSGNG(void *p) {
      delete ((::TSGNG*)p);
   }
   static void deleteArray_TSGNG(void *p) {
      delete [] ((::TSGNG*)p);
   }
   static void destruct_TSGNG(void *p) {
      typedef ::TSGNG current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TSGNG(TBuffer &buf, void *obj) {
      ((::TSGNG*)obj)->::TSGNG::Streamer(buf);
   }
} // end of namespace ROOT for class ::TSGNG

//______________________________________________________________________________
void TMLP::Streamer(TBuffer &R__b)
{
   // Stream an object of class TMLP.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TXMLP::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, TMLP::IsA());
   } else {
      R__c = R__b.WriteVersion(TMLP::IsA(), kTRUE);
      TXMLP::Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TMLP(void *p) {
      return  p ? new(p) ::TMLP : new ::TMLP;
   }
   static void *newArray_TMLP(Long_t nElements, void *p) {
      return p ? new(p) ::TMLP[nElements] : new ::TMLP[nElements];
   }
   // Wrapper around operator delete
   static void delete_TMLP(void *p) {
      delete ((::TMLP*)p);
   }
   static void deleteArray_TMLP(void *p) {
      delete [] ((::TMLP*)p);
   }
   static void destruct_TMLP(void *p) {
      typedef ::TMLP current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TMLP(TBuffer &buf, void *obj) {
      ((::TMLP*)obj)->::TMLP::Streamer(buf);
   }
} // end of namespace ROOT for class ::TMLP

//______________________________________________________________________________
void TGCS::Streamer(TBuffer &R__b)
{
   // Stream an object of class TGCS.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      VUnsupervisedNet::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, TGCS::IsA());
   } else {
      R__c = R__b.WriteVersion(TGCS::IsA(), kTRUE);
      VUnsupervisedNet::Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TGCS(void *p) {
      return  p ? new(p) ::TGCS : new ::TGCS;
   }
   static void *newArray_TGCS(Long_t nElements, void *p) {
      return p ? new(p) ::TGCS[nElements] : new ::TGCS[nElements];
   }
   // Wrapper around operator delete
   static void delete_TGCS(void *p) {
      delete ((::TGCS*)p);
   }
   static void deleteArray_TGCS(void *p) {
      delete [] ((::TGCS*)p);
   }
   static void destruct_TGCS(void *p) {
      typedef ::TGCS current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TGCS(TBuffer &buf, void *obj) {
      ((::TGCS*)obj)->::TGCS::Streamer(buf);
   }
} // end of namespace ROOT for class ::TGCS

//______________________________________________________________________________
void TGNG::Streamer(TBuffer &R__b)
{
   // Stream an object of class TGNG.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      VUnsupervisedNet::Streamer(R__b);
      R__b >> fMinDistSquare1;
      R__b >> fMinDistSquare2;
      R__b.CheckByteCount(R__s, R__c, TGNG::IsA());
   } else {
      R__c = R__b.WriteVersion(TGNG::IsA(), kTRUE);
      VUnsupervisedNet::Streamer(R__b);
      R__b << fMinDistSquare1;
      R__b << fMinDistSquare2;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TGNG(void *p) {
      return  p ? new(p) ::TGNG : new ::TGNG;
   }
   static void *newArray_TGNG(Long_t nElements, void *p) {
      return p ? new(p) ::TGNG[nElements] : new ::TGNG[nElements];
   }
   // Wrapper around operator delete
   static void delete_TGNG(void *p) {
      delete ((::TGNG*)p);
   }
   static void deleteArray_TGNG(void *p) {
      delete [] ((::TGNG*)p);
   }
   static void destruct_TGNG(void *p) {
      typedef ::TGNG current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TGNG(TBuffer &buf, void *obj) {
      ((::TGNG*)obj)->::TGNG::Streamer(buf);
   }
} // end of namespace ROOT for class ::TGNG

//______________________________________________________________________________
void TLVQ::Streamer(TBuffer &R__b)
{
   // Stream an object of class TLVQ.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      VUnsupervisedNet::Streamer(R__b);
      R__b.CheckByteCount(R__s, R__c, TLVQ::IsA());
   } else {
      R__c = R__b.WriteVersion(TLVQ::IsA(), kTRUE);
      VUnsupervisedNet::Streamer(R__b);
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TLVQ(void *p) {
      return  p ? new(p) ::TLVQ : new ::TLVQ;
   }
   static void *newArray_TLVQ(Long_t nElements, void *p) {
      return p ? new(p) ::TLVQ[nElements] : new ::TLVQ[nElements];
   }
   // Wrapper around operator delete
   static void delete_TLVQ(void *p) {
      delete ((::TLVQ*)p);
   }
   static void deleteArray_TLVQ(void *p) {
      delete [] ((::TLVQ*)p);
   }
   static void destruct_TLVQ(void *p) {
      typedef ::TLVQ current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TLVQ(TBuffer &buf, void *obj) {
      ((::TLVQ*)obj)->::TLVQ::Streamer(buf);
   }
} // end of namespace ROOT for class ::TLVQ

namespace ROOT {
   // Wrappers around operator new
   static void *new_TDataServe(void *p) {
      return  p ? new(p) ::TDataServe : new ::TDataServe;
   }
   static void *newArray_TDataServe(Long_t nElements, void *p) {
      return p ? new(p) ::TDataServe[nElements] : new ::TDataServe[nElements];
   }
   // Wrapper around operator delete
   static void delete_TDataServe(void *p) {
      delete ((::TDataServe*)p);
   }
   static void deleteArray_TDataServe(void *p) {
      delete [] ((::TDataServe*)p);
   }
   static void destruct_TDataServe(void *p) {
      typedef ::TDataServe current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TDataServe(TBuffer &buf, void *obj) {
      ((::TDataServe*)obj)->::TDataServe::Streamer(buf);
   }
} // end of namespace ROOT for class ::TDataServe

//______________________________________________________________________________
void TNNK::Streamer(TBuffer &R__b)
{
   // Stream an object of class TNNK.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      VSupervisedNet::Streamer(R__b);
      R__b >> fKernel;
      R__b.CheckByteCount(R__s, R__c, TNNK::IsA());
   } else {
      R__c = R__b.WriteVersion(TNNK::IsA(), kTRUE);
      VSupervisedNet::Streamer(R__b);
      R__b << fKernel;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TNNK(void *p) {
      return  p ? new(p) ::TNNK : new ::TNNK;
   }
   static void *newArray_TNNK(Long_t nElements, void *p) {
      return p ? new(p) ::TNNK[nElements] : new ::TNNK[nElements];
   }
   // Wrapper around operator delete
   static void delete_TNNK(void *p) {
      delete ((::TNNK*)p);
   }
   static void deleteArray_TNNK(void *p) {
      delete [] ((::TNNK*)p);
   }
   static void destruct_TNNK(void *p) {
      typedef ::TNNK current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TNNK(TBuffer &buf, void *obj) {
      ((::TNNK*)obj)->::TNNK::Streamer(buf);
   }
} // end of namespace ROOT for class ::TNNK

//______________________________________________________________________________
void TNNKernel::Streamer(TBuffer &R__b)
{
   // Stream an object of class TNNKernel.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TNamed::Streamer(R__b);
      R__b >> fNHiddL;
      R__b >> fNTrainEvents;
      R__b >> fNValidEvents;
      R__b >> fLearnParam;
      R__b >> fLowerInitWeight;
      R__b >> fUpperInitWeight;
      R__b >> fNTrainCycles;
      R__b >> fUseBiases;
      fRandom.Streamer(R__b);
      R__b >> fNWeights;
      R__b >> fMu;
      R__b >> fFlatSE;
      R__b.CheckByteCount(R__s, R__c, TNNKernel::IsA());
   } else {
      R__c = R__b.WriteVersion(TNNKernel::IsA(), kTRUE);
      TNamed::Streamer(R__b);
      R__b << fNHiddL;
      R__b << fNTrainEvents;
      R__b << fNValidEvents;
      R__b << fLearnParam;
      R__b << fLowerInitWeight;
      R__b << fUpperInitWeight;
      R__b << fNTrainCycles;
      R__b << fUseBiases;
      fRandom.Streamer(R__b);
      R__b << fNWeights;
      R__b << fMu;
      R__b << fFlatSE;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TNNKernel(void *p) {
      return  p ? new(p) ::TNNKernel : new ::TNNKernel;
   }
   static void *newArray_TNNKernel(Long_t nElements, void *p) {
      return p ? new(p) ::TNNKernel[nElements] : new ::TNNKernel[nElements];
   }
   // Wrapper around operator delete
   static void delete_TNNKernel(void *p) {
      delete ((::TNNKernel*)p);
   }
   static void deleteArray_TNNKernel(void *p) {
      delete [] ((::TNNKernel*)p);
   }
   static void destruct_TNNKernel(void *p) {
      typedef ::TNNKernel current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TNNKernel(TBuffer &buf, void *obj) {
      ((::TNNKernel*)obj)->::TNNKernel::Streamer(buf);
   }
} // end of namespace ROOT for class ::TNNKernel

//______________________________________________________________________________
void TGNGTracker::Streamer(TBuffer &R__b)
{
   // Stream an object of class TGNGTracker.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      VUnsupervisedNet::Streamer(R__b);
      R__b >> fMinDistSquare1;
      R__b >> fMinDistSquare2;
      R__b.CheckByteCount(R__s, R__c, TGNGTracker::IsA());
   } else {
      R__c = R__b.WriteVersion(TGNGTracker::IsA(), kTRUE);
      VUnsupervisedNet::Streamer(R__b);
      R__b << fMinDistSquare1;
      R__b << fMinDistSquare2;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TGNGTracker(void *p) {
      return  p ? new(p) ::TGNGTracker : new ::TGNGTracker;
   }
   static void *newArray_TGNGTracker(Long_t nElements, void *p) {
      return p ? new(p) ::TGNGTracker[nElements] : new ::TGNGTracker[nElements];
   }
   // Wrapper around operator delete
   static void delete_TGNGTracker(void *p) {
      delete ((::TGNGTracker*)p);
   }
   static void deleteArray_TGNGTracker(void *p) {
      delete [] ((::TGNGTracker*)p);
   }
   static void destruct_TGNGTracker(void *p) {
      typedef ::TGNGTracker current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TGNGTracker(TBuffer &buf, void *obj) {
      ((::TGNGTracker*)obj)->::TGNGTracker::Streamer(buf);
   }
} // end of namespace ROOT for class ::TGNGTracker

//______________________________________________________________________________
void TRadon::Streamer(TBuffer &R__b)
{
   // Stream an object of class TRadon.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> sigma;
      R__b >> threshold;
      R__b >> nt1;
      R__b >> nt2;
      Hlist.Streamer(R__b);
      {
         vector<TVector3> &R__stl =  hits;
         R__stl.clear();
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            TVector3 R__t;
            R__t.Streamer(R__b);
            R__stl.push_back(R__t);
         }
      }
      {
         vector<RADON> &R__stl =  rt;
         R__stl.clear();
         TClass *R__tcl1 = TBuffer::GetClass(typeid(RADON));
         if (R__tcl1==0) {
            Error("rt streamer","Missing the TClass object for RADON!");
            return;
         }
         int R__i, R__n;
         R__b >> R__n;
         R__stl.reserve(R__n);
         for (R__i = 0; R__i < R__n; R__i++) {
            RADON R__t;
            R__b.StreamObject(&R__t,R__tcl1);
            R__stl.push_back(R__t);
         }
      }
      R__b.CheckByteCount(R__s, R__c, TRadon::IsA());
   } else {
      R__c = R__b.WriteVersion(TRadon::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << sigma;
      R__b << threshold;
      R__b << nt1;
      R__b << nt2;
      Hlist.Streamer(R__b);
      {
         vector<TVector3> &R__stl =  hits;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
            vector<TVector3>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            ((TVector3&)(*R__k)).Streamer(R__b);
            }
         }
      }
      {
         vector<RADON> &R__stl =  rt;
         int R__n=int(R__stl.size());
         R__b << R__n;
         if(R__n) {
         TClass *R__tcl1 = TBuffer::GetClass(typeid(RADON));
         if (R__tcl1==0) {
            Error("rt streamer","Missing the TClass object for RADON!");
            return;
         }
            vector<RADON>::iterator R__k;
            for (R__k = R__stl.begin(); R__k != R__stl.end(); ++R__k) {
            R__b.StreamObject((RADON*)&(*R__k),R__tcl1);
            }
         }
      }
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TRadon(void *p) {
      return  p ? new(p) ::TRadon : new ::TRadon;
   }
   static void *newArray_TRadon(Long_t nElements, void *p) {
      return p ? new(p) ::TRadon[nElements] : new ::TRadon[nElements];
   }
   // Wrapper around operator delete
   static void delete_TRadon(void *p) {
      delete ((::TRadon*)p);
   }
   static void deleteArray_TRadon(void *p) {
      delete [] ((::TRadon*)p);
   }
   static void destruct_TRadon(void *p) {
      typedef ::TRadon current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_TRadon(TBuffer &buf, void *obj) {
      ((::TRadon*)obj)->::TRadon::Streamer(buf);
   }
} // end of namespace ROOT for class ::TRadon

namespace ROOT {
   static TClass *vectorlETVector3gR_Dictionary();
   static void vectorlETVector3gR_TClassManip(TClass*);
   static void *new_vectorlETVector3gR(void *p = 0);
   static void *newArray_vectorlETVector3gR(Long_t size, void *p);
   static void delete_vectorlETVector3gR(void *p);
   static void deleteArray_vectorlETVector3gR(void *p);
   static void destruct_vectorlETVector3gR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<TVector3>*)
   {
      vector<TVector3> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<TVector3>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<TVector3>", -2, "vector", 447,
                  typeid(vector<TVector3>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlETVector3gR_Dictionary, isa_proxy, 0,
                  sizeof(vector<TVector3>) );
      instance.SetNew(&new_vectorlETVector3gR);
      instance.SetNewArray(&newArray_vectorlETVector3gR);
      instance.SetDelete(&delete_vectorlETVector3gR);
      instance.SetDeleteArray(&deleteArray_vectorlETVector3gR);
      instance.SetDestructor(&destruct_vectorlETVector3gR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<TVector3> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<TVector3>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlETVector3gR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<TVector3>*)0x0)->GetClass();
      vectorlETVector3gR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlETVector3gR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlETVector3gR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TVector3> : new vector<TVector3>;
   }
   static void *newArray_vectorlETVector3gR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<TVector3>[nElements] : new vector<TVector3>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlETVector3gR(void *p) {
      delete ((vector<TVector3>*)p);
   }
   static void deleteArray_vectorlETVector3gR(void *p) {
      delete [] ((vector<TVector3>*)p);
   }
   static void destruct_vectorlETVector3gR(void *p) {
      typedef vector<TVector3> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<TVector3>

namespace ROOT {
   static TClass *vectorlERADONgR_Dictionary();
   static void vectorlERADONgR_TClassManip(TClass*);
   static void *new_vectorlERADONgR(void *p = 0);
   static void *newArray_vectorlERADONgR(Long_t size, void *p);
   static void delete_vectorlERADONgR(void *p);
   static void deleteArray_vectorlERADONgR(void *p);
   static void destruct_vectorlERADONgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<RADON>*)
   {
      vector<RADON> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<RADON>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<RADON>", -2, "vector", 447,
                  typeid(vector<RADON>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlERADONgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<RADON>) );
      instance.SetNew(&new_vectorlERADONgR);
      instance.SetNewArray(&newArray_vectorlERADONgR);
      instance.SetDelete(&delete_vectorlERADONgR);
      instance.SetDeleteArray(&deleteArray_vectorlERADONgR);
      instance.SetDestructor(&destruct_vectorlERADONgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<RADON> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<RADON>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlERADONgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<RADON>*)0x0)->GetClass();
      vectorlERADONgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlERADONgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlERADONgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RADON> : new vector<RADON>;
   }
   static void *newArray_vectorlERADONgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<RADON>[nElements] : new vector<RADON>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlERADONgR(void *p) {
      delete ((vector<RADON>*)p);
   }
   static void deleteArray_vectorlERADONgR(void *p) {
      delete [] ((vector<RADON>*)p);
   }
   static void destruct_vectorlERADONgR(void *p) {
      typedef vector<RADON> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<RADON>

namespace {
  void TriggerDictionaryInitialization_RhoNNOCint_Impl() {
    static const char* headers[] = {
"VNeuralNetPlotter.h",
"VNeuralNet.h",
"VSupervisedNet.h",
"VUnsupervisedNet.h",
"TFD.h",
"TSGCS.h",
"TSGNG.h",
"TPerceptron.h",
"TMLP.h",
"TXMLP.h",
"TNeuralNetCell.h",
"TGCS.h",
"TGNG.h",
"TLVQ.h",
"TDataServe.h",
"TNNK.h",
"TGNGTracker.h",
"TRadon.h",
0
    };
    static const char* includePaths[] = {
"..",
"/Applications/root_v6.12.06/include",
"/Users/marcel/workspace/rhonno/RhoNNO/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "RhoNNOCint dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
class __attribute__((annotate(R"ATTRDUMP(Base class of all network plotters)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$VNeuralNetPlotter.h")))  VNeuralNetPlotter;
class __attribute__((annotate(R"ATTRDUMP(Base class of all network plotters)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$VNeuralNetPlotter.h")))  TSimpleNeuralNetPlotter;
class __attribute__((annotate(R"ATTRDUMP(Parameters for all supervised networks)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$VNeuralNet.h")))  TNeuralNetParameters;
class __attribute__((annotate(R"ATTRDUMP(Base class of all networks)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$VNeuralNet.h")))  VNeuralNet;
class __attribute__((annotate(R"ATTRDUMP(Supervised training)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$VSupervisedNet.h")))  VSupervisedNet;
class __attribute__((annotate(R"ATTRDUMP(Parameters for unsupervised networks)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$RhoNNO/TNeuralNetCell.h")))  __attribute__((annotate("$clingAutoload$VUnsupervisedNet.h")))  TNeuralNetCellParameters;
class __attribute__((annotate(R"ATTRDUMP(Cell for unsupervised networks)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$RhoNNO/TNeuralNetCell.h")))  __attribute__((annotate("$clingAutoload$VUnsupervisedNet.h")))  TNeuralNetCell;
class __attribute__((annotate(R"ATTRDUMP(Unsupervised training)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$VUnsupervisedNet.h")))  VUnsupervisedNet;
class __attribute__((annotate(R"ATTRDUMP(Multilayer Perceptron)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$RhoNNO/TPerceptron.h")))  __attribute__((annotate("$clingAutoload$TFD.h")))  TPerceptron;
class __attribute__((annotate(R"ATTRDUMP(Multi-Layer Perceptron (extended))ATTRDUMP"))) __attribute__((annotate("$clingAutoload$RhoNNO/TXMLP.h")))  __attribute__((annotate("$clingAutoload$TFD.h")))  TXMLP;
class __attribute__((annotate(R"ATTRDUMP(Fisher discriminant)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TFD.h")))  TFD;
class __attribute__((annotate(R"ATTRDUMP(Supervised Growing Cell Structure)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TSGCS.h")))  TSGCS;
class __attribute__((annotate(R"ATTRDUMP(Supervised Growing Neural Gas)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TSGNG.h")))  TSGNG;
class __attribute__((annotate(R"ATTRDUMP(Multi-Layer Perceptron)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TMLP.h")))  TMLP;
class __attribute__((annotate(R"ATTRDUMP(Unsupervised Growing Cellstructure)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TGCS.h")))  TGCS;
class __attribute__((annotate(R"ATTRDUMP(Unsupervised Growing Neural Gas)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TGNG.h")))  TGNG;
class __attribute__((annotate(R"ATTRDUMP(Learning Vector Quantization)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TLVQ.h")))  TLVQ;
class __attribute__((annotate(R"ATTRDUMP(Database for network training)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TDataServe.h")))  TDataServe;
class __attribute__((annotate(R"ATTRDUMP(Interface to TNNKernel)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TNNK.h")))  TNNK;
class __attribute__((annotate("$clingAutoload$TNNK.h")))  TNNKernel;
class __attribute__((annotate(R"ATTRDUMP(Unsupervised Growing Neural Gas)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TGNGTracker.h")))  TGNGTracker;
class __attribute__((annotate(R"ATTRDUMP(Fuzzy Radon transform)ATTRDUMP"))) __attribute__((annotate("$clingAutoload$TRadon.h")))  TRadon;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "RhoNNOCint dictionary payload"

#ifndef G__VECTOR_HAS_CLASS_ITERATOR
  #define G__VECTOR_HAS_CLASS_ITERATOR 1
#endif

#define _BACKWARD_BACKWARD_WARNING_H
#include "VNeuralNetPlotter.h"
#include "VNeuralNet.h"
#include "VSupervisedNet.h"
#include "VUnsupervisedNet.h"
#include "TFD.h"
#include "TSGCS.h"
#include "TSGNG.h"
#include "TPerceptron.h"
#include "TMLP.h"
#include "TXMLP.h"
#include "TNeuralNetCell.h"
#include "TGCS.h"
#include "TGNG.h"
#include "TLVQ.h"
#include "TDataServe.h"
#include "TNNK.h"
#include "TGNGTracker.h"
#include "TRadon.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[]={
"TDataServe", payloadCode, "@",
"TFD", payloadCode, "@",
"TGCS", payloadCode, "@",
"TGNG", payloadCode, "@",
"TGNGTracker", payloadCode, "@",
"TLVQ", payloadCode, "@",
"TMLP", payloadCode, "@",
"TNNK", payloadCode, "@",
"TNNKernel", payloadCode, "@",
"TNeuralNetCell", payloadCode, "@",
"TNeuralNetCellParameters", payloadCode, "@",
"TNeuralNetParameters", payloadCode, "@",
"TPerceptron", payloadCode, "@",
"TRadon", payloadCode, "@",
"TSGCS", payloadCode, "@",
"TSGNG", payloadCode, "@",
"TSimpleNeuralNetPlotter", payloadCode, "@",
"TXMLP", payloadCode, "@",
"VNeuralNet", payloadCode, "@",
"VNeuralNetPlotter", payloadCode, "@",
"VSupervisedNet", payloadCode, "@",
"VUnsupervisedNet", payloadCode, "@",
nullptr};

    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("RhoNNOCint",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_RhoNNOCint_Impl, {}, classesHeaders);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_RhoNNOCint_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_RhoNNOCint() {
  TriggerDictionaryInitialization_RhoNNOCint_Impl();
}
