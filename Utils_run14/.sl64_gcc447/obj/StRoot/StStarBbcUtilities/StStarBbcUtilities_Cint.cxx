//
// File generated by rootcint at Mon Mar 23 10:05:12 2015

// Do NOT change. Changes will be lost next time file is generated
//

#define R__DICTIONARY_FILENAME dOsl64_gcc447dIobjdIStRootdIStStarBbcUtilitiesdIStStarBbcUtilities_Cint
#include "RConfig.h" //rootcint 4834
#if !defined(R__ACCESS_IN_SYMBOL)
//Break the privacy of classes -- Disabled for the moment
#define private public
#define protected public
#endif

// Since CINT ignores the std namespace, we need to do so in this file.
namespace std {} using namespace std;
#include "StStarBbcUtilities_Cint.h"

#include "TClass.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"

// START OF SHADOWS

namespace ROOT {
   namespace Shadow {
   } // of namespace Shadow
} // of namespace ROOT
// END OF SHADOWS

namespace ROOT {
   void StStarBbcUtilities_ShowMembers(void *obj, TMemberInspector &R__insp);
   static void *new_StStarBbcUtilities(void *p = 0);
   static void *newArray_StStarBbcUtilities(Long_t size, void *p);
   static void delete_StStarBbcUtilities(void *p);
   static void deleteArray_StStarBbcUtilities(void *p);
   static void destruct_StStarBbcUtilities(void *p);
   static void streamer_StStarBbcUtilities(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::StStarBbcUtilities*)
   {
      ::StStarBbcUtilities *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::StStarBbcUtilities >(0);
      static ::ROOT::TGenericClassInfo 
         instance("StStarBbcUtilities", ::StStarBbcUtilities::Class_Version(), ".sl64_gcc447/obj/StRoot/StStarBbcUtilities/StStarBbcUtilities.h", 10,
                  typeid(::StStarBbcUtilities), DefineBehavior(ptr, ptr),
                  &::StStarBbcUtilities::Dictionary, isa_proxy, 0,
                  sizeof(::StStarBbcUtilities) );
      instance.SetNew(&new_StStarBbcUtilities);
      instance.SetNewArray(&newArray_StStarBbcUtilities);
      instance.SetDelete(&delete_StStarBbcUtilities);
      instance.SetDeleteArray(&deleteArray_StStarBbcUtilities);
      instance.SetDestructor(&destruct_StStarBbcUtilities);
      instance.SetStreamerFunc(&streamer_StStarBbcUtilities);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::StStarBbcUtilities*)
   {
      return GenerateInitInstanceLocal((::StStarBbcUtilities*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_(Init) = GenerateInitInstanceLocal((const ::StStarBbcUtilities*)0x0); R__UseDummy(_R__UNIQUE_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
TClass *StStarBbcUtilities::fgIsA = 0;  // static to hold class pointer

//______________________________________________________________________________
const char *StStarBbcUtilities::Class_Name()
{
   return "StStarBbcUtilities";
}

//______________________________________________________________________________
const char *StStarBbcUtilities::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StStarBbcUtilities*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int StStarBbcUtilities::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::StStarBbcUtilities*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
void StStarBbcUtilities::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StStarBbcUtilities*)0x0)->GetClass();
}

//______________________________________________________________________________
TClass *StStarBbcUtilities::Class()
{
   if (!fgIsA) fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::StStarBbcUtilities*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
void StStarBbcUtilities::Streamer(TBuffer &R__b)
{
   // Stream an object of class StStarBbcUtilities.

   ::Error("StStarBbcUtilities::Streamer", "version id <=0 in ClassDef, dummy Streamer() called"); if (R__b.IsReading()) { }
}

//______________________________________________________________________________
void StStarBbcUtilities::ShowMembers(TMemberInspector &R__insp)
{
      // Inspect the data members of an object of class StStarBbcUtilities.
      TClass *R__cl = ::StStarBbcUtilities::IsA();
      if (R__cl || R__insp.IsA()) { }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_StStarBbcUtilities(void *p) {
      return  p ? ::new((::ROOT::TOperatorNewHelper*)p) ::StStarBbcUtilities : new ::StStarBbcUtilities;
   }
   static void *newArray_StStarBbcUtilities(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::TOperatorNewHelper*)p) ::StStarBbcUtilities[nElements] : new ::StStarBbcUtilities[nElements];
   }
   // Wrapper around operator delete
   static void delete_StStarBbcUtilities(void *p) {
      delete ((::StStarBbcUtilities*)p);
   }
   static void deleteArray_StStarBbcUtilities(void *p) {
      delete [] ((::StStarBbcUtilities*)p);
   }
   static void destruct_StStarBbcUtilities(void *p) {
      typedef ::StStarBbcUtilities current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_StStarBbcUtilities(TBuffer &buf, void *obj) {
      ((::StStarBbcUtilities*)obj)->::StStarBbcUtilities::Streamer(buf);
   }
} // end of namespace ROOT for class ::StStarBbcUtilities

/********************************************************
* .sl64_gcc447/obj/StRoot/StStarBbcUtilities/StStarBbcUtilities_Cint.cxx
* CAUTION: DON'T CHANGE THIS FILE. THIS FILE IS AUTOMATICALLY GENERATED
*          FROM HEADER FILES LISTED IN G__setup_cpp_environmentXXX().
*          CHANGE THOSE HEADER FILES AND REGENERATE THIS FILE.
********************************************************/

#ifdef G__MEMTEST
#undef malloc
#undef free
#endif

#if defined(__GNUC__) && __GNUC__ >= 4 && ((__GNUC_MINOR__ == 2 && __GNUC_PATCHLEVEL__ >= 1) || (__GNUC_MINOR__ >= 3))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

extern "C" void G__cpp_reset_tagtableStStarBbcUtilities_Cint();

extern "C" void G__set_cpp_environmentStStarBbcUtilities_Cint() {
  G__add_compiledheader("TObject.h");
  G__add_compiledheader("TMemberInspector.h");
  G__add_compiledheader("StStarBbcUtilities.h");
  G__cpp_reset_tagtableStStarBbcUtilities_Cint();
}
#include <new>
extern "C" int G__cpp_dllrevStStarBbcUtilities_Cint() { return(30051515); }

/*********************************************************
* Member function Interface Method
*********************************************************/

/* StStarBbcUtilities */
static int G__StStarBbcUtilities_Cint_168_0_1(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   StStarBbcUtilities* p = NULL;
   char* gvp = (char*) G__getgvp();
   int n = G__getaryconstruct();
   if (n) {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new StStarBbcUtilities[n];
     } else {
       p = new((void*) gvp) StStarBbcUtilities[n];
     }
   } else {
     if ((gvp == (char*)G__PVOID) || (gvp == 0)) {
       p = new StStarBbcUtilities;
     } else {
       p = new((void*) gvp) StStarBbcUtilities;
     }
   }
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__StStarBbcUtilities_CintLN_StStarBbcUtilities));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StStarBbcUtilities_Cint_168_0_2(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 114, (long) ((const StStarBbcUtilities*) G__getstructoffset())->GetAdcCut(*((const TString*) G__int(libp->para[0]))));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StStarBbcUtilities_Cint_168_0_3(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 102, (double) ((const StStarBbcUtilities*) G__getstructoffset())->GetPhi((const Int_t) G__int(libp->para[0]), (const Int_t) G__int(libp->para[1])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StStarBbcUtilities_Cint_168_0_4(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letdouble(result7, 102, (double) ((const StStarBbcUtilities*) G__getstructoffset())->GetXY((const Int_t) G__int(libp->para[0]), (const Int_t) G__int(libp->para[1])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StStarBbcUtilities_Cint_168_0_5(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((const StStarBbcUtilities*) G__getstructoffset())->GetRandomXY((const Int_t) G__int(libp->para[0]), *(Float_t*) G__Floatref(&libp->para[1])
, *(Float_t*) G__Floatref(&libp->para[2]));
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StStarBbcUtilities_Cint_168_0_6(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 104, (long) ((const StStarBbcUtilities*) G__getstructoffset())->EastWest((const UInt_t) G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StStarBbcUtilities_Cint_168_0_7(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 104, (long) ((const StStarBbcUtilities*) G__getstructoffset())->PmtId((const UInt_t) G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StStarBbcUtilities_Cint_168_0_8(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 104, (long) ((const StStarBbcUtilities*) G__getstructoffset())->InnerOuter((const UInt_t) G__int(libp->para[0])));
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StStarBbcUtilities_Cint_168_0_9(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) StStarBbcUtilities::Class());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StStarBbcUtilities_Cint_168_0_10(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) StStarBbcUtilities::Class_Name());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StStarBbcUtilities_Cint_168_0_11(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 115, (long) StStarBbcUtilities::Class_Version());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StStarBbcUtilities_Cint_168_0_12(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      StStarBbcUtilities::Dictionary();
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StStarBbcUtilities_Cint_168_0_13(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 85, (long) ((const StStarBbcUtilities*) G__getstructoffset())->IsA());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StStarBbcUtilities_Cint_168_0_14(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((StStarBbcUtilities*) G__getstructoffset())->ShowMembers(*(TMemberInspector*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StStarBbcUtilities_Cint_168_0_15(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((StStarBbcUtilities*) G__getstructoffset())->Streamer(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StStarBbcUtilities_Cint_168_0_16(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      ((StStarBbcUtilities*) G__getstructoffset())->StreamerNVirtual(*(TBuffer*) libp->para[0].ref);
      G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StStarBbcUtilities_Cint_168_0_17(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) StStarBbcUtilities::DeclFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StStarBbcUtilities_Cint_168_0_18(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) StStarBbcUtilities::ImplFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StStarBbcUtilities_Cint_168_0_19(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 67, (long) StStarBbcUtilities::ImplFileName());
   return(1 || funcname || hash || result7 || libp) ;
}

static int G__StStarBbcUtilities_Cint_168_0_20(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
      G__letint(result7, 105, (long) StStarBbcUtilities::DeclFileLine());
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic copy constructor
static int G__StStarBbcUtilities_Cint_168_0_21(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)

{
   StStarBbcUtilities* p;
   void* tmp = (void*) G__int(libp->para[0]);
   p = new StStarBbcUtilities(*(StStarBbcUtilities*) tmp);
   result7->obj.i = (long) p;
   result7->ref = (long) p;
   G__set_tagnum(result7,G__get_linked_tagnum(&G__StStarBbcUtilities_CintLN_StStarBbcUtilities));
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic destructor
typedef StStarBbcUtilities G__TStStarBbcUtilities;
static int G__StStarBbcUtilities_Cint_168_0_22(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   char* gvp = (char*) G__getgvp();
   long soff = G__getstructoffset();
   int n = G__getaryconstruct();
   //
   //has_a_delete: 0
   //has_own_delete1arg: 0
   //has_own_delete2arg: 0
   //
   if (!soff) {
     return(1);
   }
   if (n) {
     if (gvp == (char*)G__PVOID) {
       delete[] (StStarBbcUtilities*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       for (int i = n - 1; i >= 0; --i) {
         ((StStarBbcUtilities*) (soff+(sizeof(StStarBbcUtilities)*i)))->~G__TStStarBbcUtilities();
       }
       G__setgvp((long)gvp);
     }
   } else {
     if (gvp == (char*)G__PVOID) {
       delete (StStarBbcUtilities*) soff;
     } else {
       G__setgvp((long) G__PVOID);
       ((StStarBbcUtilities*) (soff))->~G__TStStarBbcUtilities();
       G__setgvp((long)gvp);
     }
   }
   G__setnull(result7);
   return(1 || funcname || hash || result7 || libp) ;
}

// automatic assignment operator
static int G__StStarBbcUtilities_Cint_168_0_23(G__value* result7, G__CONST char* funcname, struct G__param* libp, int hash)
{
   StStarBbcUtilities* dest = (StStarBbcUtilities*) G__getstructoffset();
   *dest = *(StStarBbcUtilities*) libp->para[0].ref;
   const StStarBbcUtilities& obj = *dest;
   result7->ref = (long) (&obj);
   result7->obj.i = (long) (&obj);
   return(1 || funcname || hash || result7 || libp) ;
}


/* Setting up global function */

/*********************************************************
* Member function Stub
*********************************************************/

/* StStarBbcUtilities */

/*********************************************************
* Global function Stub
*********************************************************/

/*********************************************************
* Get size of pointer to member function
*********************************************************/
class G__Sizep2memfuncStStarBbcUtilities_Cint {
 public:
  G__Sizep2memfuncStStarBbcUtilities_Cint(): p(&G__Sizep2memfuncStStarBbcUtilities_Cint::sizep2memfunc) {}
    size_t sizep2memfunc() { return(sizeof(p)); }
  private:
    size_t (G__Sizep2memfuncStStarBbcUtilities_Cint::*p)();
};

size_t G__get_sizep2memfuncStStarBbcUtilities_Cint()
{
  G__Sizep2memfuncStStarBbcUtilities_Cint a;
  G__setsizep2memfunc((int)a.sizep2memfunc());
  return((size_t)a.sizep2memfunc());
}


/*********************************************************
* virtual base class offset calculation interface
*********************************************************/

   /* Setting up class inheritance */

/*********************************************************
* Inheritance information setup/
*********************************************************/
extern "C" void G__cpp_setup_inheritanceStStarBbcUtilities_Cint() {

   /* Setting up class inheritance */
}

/*********************************************************
* typedef information setup/
*********************************************************/
extern "C" void G__cpp_setup_typetableStStarBbcUtilities_Cint() {

   /* Setting up typedef entry */
   G__search_typename2("UShort_t",114,-1,0,-1);
   G__setnewtype(-1,"Unsigned Short integer 2 bytes (unsigned short)",0);
   G__search_typename2("Int_t",105,-1,0,-1);
   G__setnewtype(-1,"Signed integer 4 bytes (int)",0);
   G__search_typename2("UInt_t",104,-1,0,-1);
   G__setnewtype(-1,"Unsigned integer 4 bytes (unsigned int)",0);
   G__search_typename2("Float_t",102,-1,0,-1);
   G__setnewtype(-1,"Float 4 bytes (float)",0);
   G__search_typename2("Version_t",115,-1,0,-1);
   G__setnewtype(-1,"Class version identifier (short)",0);
   G__search_typename2("vector<ROOT::TSchemaHelper>",117,G__get_linked_tagnum(&G__StStarBbcUtilities_CintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__StStarBbcUtilities_CintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__StStarBbcUtilities_CintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__StStarBbcUtilities_CintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__StStarBbcUtilities_CintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("vector<TVirtualArray*>",117,G__get_linked_tagnum(&G__StStarBbcUtilities_CintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR),0,-1);
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<const_iterator>",117,G__get_linked_tagnum(&G__StStarBbcUtilities_CintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__StStarBbcUtilities_CintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
   G__search_typename2("reverse_iterator<iterator>",117,G__get_linked_tagnum(&G__StStarBbcUtilities_CintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR),0,G__get_linked_tagnum(&G__StStarBbcUtilities_CintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR));
   G__setnewtype(-1,NULL,0);
}

/*********************************************************
* Data Member information setup/
*********************************************************/

   /* Setting up class,struct,union tag member variable */

   /* StStarBbcUtilities */
static void G__setup_memvarStStarBbcUtilities(void) {
   G__tag_memvar_setup(G__get_linked_tagnum(&G__StStarBbcUtilities_CintLN_StStarBbcUtilities));
   { StStarBbcUtilities *p; p=(StStarBbcUtilities*)0x1000; if (p) { }
   G__memvar_setup((void*)0,108,0,0,-1,-1,-1,4,"G__virtualinfo=",0,(char*)NULL);
   G__memvar_setup((void*)0,85,0,0,G__get_linked_tagnum(&G__StStarBbcUtilities_CintLN_TClass),-1,-2,4,"fgIsA=",0,(char*)NULL);
   }
   G__tag_memvar_reset();
}

extern "C" void G__cpp_setup_memvarStStarBbcUtilities_Cint() {
}
/***********************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
************************************************************
***********************************************************/

/*********************************************************
* Member function information setup for each class
*********************************************************/
static void G__setup_memfuncStStarBbcUtilities(void) {
   /* StStarBbcUtilities */
   G__tag_memfunc_setup(G__get_linked_tagnum(&G__StStarBbcUtilities_CintLN_StStarBbcUtilities));
   G__memfunc_setup("StStarBbcUtilities",1828,G__StStarBbcUtilities_Cint_168_0_1, 105, G__get_linked_tagnum(&G__StStarBbcUtilities_CintLN_StStarBbcUtilities), -1, 0, 0, 1, 1, 0, "", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetAdcCut",852,G__StStarBbcUtilities_Cint_168_0_2, 114, -1, G__defined_typename("UShort_t"), 0, 1, 1, 1, 8, "u 'TString' - 10 - energy", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetPhi",577,G__StStarBbcUtilities_Cint_168_0_3, 102, -1, G__defined_typename("Float_t"), 0, 2, 1, 1, 8, 
"i - 'Int_t' 10 - eastWest i - 'Int_t' 10 - tileId", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetXY",465,G__StStarBbcUtilities_Cint_168_0_4, 102, -1, G__defined_typename("Float_t"), 0, 2, 1, 1, 8, 
"i - 'Int_t' 10 - tileId i - 'Int_t' 10 - idXY", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("GetRandomXY",1074,G__StStarBbcUtilities_Cint_168_0_5, 121, -1, -1, 0, 3, 1, 1, 8, 
"i - 'Int_t' 10 - tileId f - 'Float_t' 1 - x_pos_random "
"f - 'Float_t' 1 - y_pos_random", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("EastWest",816,G__StStarBbcUtilities_Cint_168_0_6, 104, -1, G__defined_typename("UInt_t"), 0, 1, 1, 1, 8, "h - 'UInt_t' 10 - ipmt_all", "/ East:0-23, West:24-47", (void*) NULL, 0);
   G__memfunc_setup("PmtId",478,G__StStarBbcUtilities_Cint_168_0_7, 104, -1, G__defined_typename("UInt_t"), 0, 1, 1, 1, 8, "h - 'UInt_t' 10 - ipmt_all", "/ PMT id for each arm: 0-23", (void*) NULL, 0);
   G__memfunc_setup("InnerOuter",1035,G__StStarBbcUtilities_Cint_168_0_8, 104, -1, G__defined_typename("UInt_t"), 0, 1, 1, 1, 8, "h - 'UInt_t' 10 - ipmt_each", "/ Innter:0-15, Outer:16-23", (void*) NULL, 0);
   G__memfunc_setup("Class",502,G__StStarBbcUtilities_Cint_168_0_9, 85, G__get_linked_tagnum(&G__StStarBbcUtilities_CintLN_TClass), -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (TClass* (*)())(&StStarBbcUtilities::Class) ), 0);
   G__memfunc_setup("Class_Name",982,G__StStarBbcUtilities_Cint_168_0_10, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&StStarBbcUtilities::Class_Name) ), 0);
   G__memfunc_setup("Class_Version",1339,G__StStarBbcUtilities_Cint_168_0_11, 115, -1, G__defined_typename("Version_t"), 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (Version_t (*)())(&StStarBbcUtilities::Class_Version) ), 0);
   G__memfunc_setup("Dictionary",1046,G__StStarBbcUtilities_Cint_168_0_12, 121, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (void (*)())(&StStarBbcUtilities::Dictionary) ), 0);
   G__memfunc_setup("IsA",253,G__StStarBbcUtilities_Cint_168_0_13, 85, G__get_linked_tagnum(&G__StStarBbcUtilities_CintLN_TClass), -1, 0, 0, 1, 1, 8, "", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("ShowMembers",1132,G__StStarBbcUtilities_Cint_168_0_14, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TMemberInspector' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("Streamer",835,G__StStarBbcUtilities_Cint_168_0_15, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - -", (char*)NULL, (void*) NULL, 1);
   G__memfunc_setup("StreamerNVirtual",1656,G__StStarBbcUtilities_Cint_168_0_16, 121, -1, -1, 0, 1, 1, 1, 0, "u 'TBuffer' - 1 - ClassDef_StreamerNVirtual_b", (char*)NULL, (void*) NULL, 0);
   G__memfunc_setup("DeclFileName",1145,G__StStarBbcUtilities_Cint_168_0_17, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&StStarBbcUtilities::DeclFileName) ), 0);
   G__memfunc_setup("ImplFileLine",1178,G__StStarBbcUtilities_Cint_168_0_18, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&StStarBbcUtilities::ImplFileLine) ), 0);
   G__memfunc_setup("ImplFileName",1171,G__StStarBbcUtilities_Cint_168_0_19, 67, -1, -1, 0, 0, 3, 1, 1, "", (char*)NULL, (void*) G__func2void( (const char* (*)())(&StStarBbcUtilities::ImplFileName) ), 0);
   G__memfunc_setup("DeclFileLine",1152,G__StStarBbcUtilities_Cint_168_0_20, 105, -1, -1, 0, 0, 3, 1, 0, "", (char*)NULL, (void*) G__func2void( (int (*)())(&StStarBbcUtilities::DeclFileLine) ), 0);
   // automatic copy constructor
   G__memfunc_setup("StStarBbcUtilities", 1828, G__StStarBbcUtilities_Cint_168_0_21, (int) ('i'), G__get_linked_tagnum(&G__StStarBbcUtilities_CintLN_StStarBbcUtilities), -1, 0, 1, 1, 1, 0, "u 'StStarBbcUtilities' - 11 - -", (char*) NULL, (void*) NULL, 0);
   // automatic destructor
   G__memfunc_setup("~StStarBbcUtilities", 1954, G__StStarBbcUtilities_Cint_168_0_22, (int) ('y'), -1, -1, 0, 0, 1, 1, 0, "", (char*) NULL, (void*) NULL, 1);
   // automatic assignment operator
   G__memfunc_setup("operator=", 937, G__StStarBbcUtilities_Cint_168_0_23, (int) ('u'), G__get_linked_tagnum(&G__StStarBbcUtilities_CintLN_StStarBbcUtilities), -1, 1, 1, 1, 1, 0, "u 'StStarBbcUtilities' - 11 - -", (char*) NULL, (void*) NULL, 0);
   G__tag_memfunc_reset();
}


/*********************************************************
* Member function information setup
*********************************************************/
extern "C" void G__cpp_setup_memfuncStStarBbcUtilities_Cint() {
}

/*********************************************************
* Global variable information setup for each class
*********************************************************/
static void G__cpp_setup_global0() {

   /* Setting up global variables */
   G__resetplocal();

}

static void G__cpp_setup_global1() {

   G__resetglobalenv();
}
extern "C" void G__cpp_setup_globalStStarBbcUtilities_Cint() {
  G__cpp_setup_global0();
  G__cpp_setup_global1();
}

/*********************************************************
* Global function information setup for each class
*********************************************************/
static void G__cpp_setup_func0() {
   G__lastifuncposition();

}

static void G__cpp_setup_func1() {
}

static void G__cpp_setup_func2() {
}

static void G__cpp_setup_func3() {
}

static void G__cpp_setup_func4() {
}

static void G__cpp_setup_func5() {
}

static void G__cpp_setup_func6() {
}

static void G__cpp_setup_func7() {
}

static void G__cpp_setup_func8() {
}

static void G__cpp_setup_func9() {
}

static void G__cpp_setup_func10() {
}

static void G__cpp_setup_func11() {
}

static void G__cpp_setup_func12() {

   G__resetifuncposition();
}

extern "C" void G__cpp_setup_funcStStarBbcUtilities_Cint() {
  G__cpp_setup_func0();
  G__cpp_setup_func1();
  G__cpp_setup_func2();
  G__cpp_setup_func3();
  G__cpp_setup_func4();
  G__cpp_setup_func5();
  G__cpp_setup_func6();
  G__cpp_setup_func7();
  G__cpp_setup_func8();
  G__cpp_setup_func9();
  G__cpp_setup_func10();
  G__cpp_setup_func11();
  G__cpp_setup_func12();
}

/*********************************************************
* Class,struct,union,enum tag information setup
*********************************************************/
/* Setup class/struct taginfo */
G__linked_taginfo G__StStarBbcUtilities_CintLN_TClass = { "TClass" , 99 , -1 };
G__linked_taginfo G__StStarBbcUtilities_CintLN_TBuffer = { "TBuffer" , 99 , -1 };
G__linked_taginfo G__StStarBbcUtilities_CintLN_TMemberInspector = { "TMemberInspector" , 99 , -1 };
G__linked_taginfo G__StStarBbcUtilities_CintLN_TString = { "TString" , 99 , -1 };
G__linked_taginfo G__StStarBbcUtilities_CintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR = { "vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >" , 99 , -1 };
G__linked_taginfo G__StStarBbcUtilities_CintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR = { "reverse_iterator<vector<ROOT::TSchemaHelper,allocator<ROOT::TSchemaHelper> >::iterator>" , 99 , -1 };
G__linked_taginfo G__StStarBbcUtilities_CintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR = { "vector<TVirtualArray*,allocator<TVirtualArray*> >" , 99 , -1 };
G__linked_taginfo G__StStarBbcUtilities_CintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR = { "reverse_iterator<vector<TVirtualArray*,allocator<TVirtualArray*> >::iterator>" , 99 , -1 };
G__linked_taginfo G__StStarBbcUtilities_CintLN_StStarBbcUtilities = { "StStarBbcUtilities" , 99 , -1 };

/* Reset class/struct taginfo */
extern "C" void G__cpp_reset_tagtableStStarBbcUtilities_Cint() {
  G__StStarBbcUtilities_CintLN_TClass.tagnum = -1 ;
  G__StStarBbcUtilities_CintLN_TBuffer.tagnum = -1 ;
  G__StStarBbcUtilities_CintLN_TMemberInspector.tagnum = -1 ;
  G__StStarBbcUtilities_CintLN_TString.tagnum = -1 ;
  G__StStarBbcUtilities_CintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR.tagnum = -1 ;
  G__StStarBbcUtilities_CintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__StStarBbcUtilities_CintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR.tagnum = -1 ;
  G__StStarBbcUtilities_CintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR.tagnum = -1 ;
  G__StStarBbcUtilities_CintLN_StStarBbcUtilities.tagnum = -1 ;
}


extern "C" void G__cpp_setup_tagtableStStarBbcUtilities_Cint() {

   /* Setting up class,struct,union tag entry */
   G__get_linked_tagnum_fwd(&G__StStarBbcUtilities_CintLN_TClass);
   G__get_linked_tagnum_fwd(&G__StStarBbcUtilities_CintLN_TBuffer);
   G__get_linked_tagnum_fwd(&G__StStarBbcUtilities_CintLN_TMemberInspector);
   G__get_linked_tagnum_fwd(&G__StStarBbcUtilities_CintLN_TString);
   G__get_linked_tagnum_fwd(&G__StStarBbcUtilities_CintLN_vectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgR);
   G__get_linked_tagnum_fwd(&G__StStarBbcUtilities_CintLN_reverse_iteratorlEvectorlEROOTcLcLTSchemaHelpercOallocatorlEROOTcLcLTSchemaHelpergRsPgRcLcLiteratorgR);
   G__get_linked_tagnum_fwd(&G__StStarBbcUtilities_CintLN_vectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgR);
   G__get_linked_tagnum_fwd(&G__StStarBbcUtilities_CintLN_reverse_iteratorlEvectorlETVirtualArraymUcOallocatorlETVirtualArraymUgRsPgRcLcLiteratorgR);
   G__tagtable_setup(G__get_linked_tagnum_fwd(&G__StStarBbcUtilities_CintLN_StStarBbcUtilities),sizeof(StStarBbcUtilities),-1,1280,(char*)NULL,G__setup_memvarStStarBbcUtilities,G__setup_memfuncStStarBbcUtilities);
}
extern "C" void G__cpp_setupStStarBbcUtilities_Cint(void) {
  G__check_setup_version(30051515,"G__cpp_setupStStarBbcUtilities_Cint()");
  G__set_cpp_environmentStStarBbcUtilities_Cint();
  G__cpp_setup_tagtableStStarBbcUtilities_Cint();

  G__cpp_setup_inheritanceStStarBbcUtilities_Cint();

  G__cpp_setup_typetableStStarBbcUtilities_Cint();

  G__cpp_setup_memvarStStarBbcUtilities_Cint();

  G__cpp_setup_memfuncStStarBbcUtilities_Cint();
  G__cpp_setup_globalStStarBbcUtilities_Cint();
  G__cpp_setup_funcStStarBbcUtilities_Cint();

   if(0==G__getsizep2memfunc()) G__get_sizep2memfuncStStarBbcUtilities_Cint();
  return;
}
class G__cpp_setup_initStStarBbcUtilities_Cint {
  public:
    G__cpp_setup_initStStarBbcUtilities_Cint() { G__add_setup_func("StStarBbcUtilities_Cint",(G__incsetup)(&G__cpp_setupStStarBbcUtilities_Cint)); G__call_setup_funcs(); }
   ~G__cpp_setup_initStStarBbcUtilities_Cint() { G__remove_setup_func("StStarBbcUtilities_Cint"); }
};
G__cpp_setup_initStStarBbcUtilities_Cint G__cpp_setup_initializerStStarBbcUtilities_Cint;

