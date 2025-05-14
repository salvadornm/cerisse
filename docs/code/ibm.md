# IBM Solver

This page explains how the [IBM method](../theory/ibmeb.md) is implemented in the code. **Cerisse** uses two boolean arrays to modify the numerical solvers

```cpp
        ibMarkers(i, j, k, 0) ;  
        ibMarkers(i, j, k, 1) ;    
```

The first boolean `ibMarkers(i, j, k, 0)` indicates if the cell is solid or not, if **true**, then is a solid cell or internal (**SP**). The second boolean `ibMarkers(i, j, k, 1)` if **true**, shows that the cell is a _ghost point_ (**GP**) and therefore requires especial treatment. A ghost point is a solid point.

## Interpolation

## Boundary Condition

### Options

### User-specific boundaries

Is possible to create a fully specific wall bc in IBM using templates in `prob.h`\
Following the exampe of manual sources or user-specifc EB boundaries

1. Declare a new class : For example `ibm_user_t`

```cpp
template < typename param, typename cls_t> class ibm_user_t;
```

The only parameter needed would be the class `cls_t`, the rest is user defined\


2. Use this class as bc as usual

```cpp
typedef ibm_user_t<ProbParm, ProbClosures> TypeWall;
typedef eib_t<TypeWall, ibm_param_t ProbClosures> ProbIB;
```

&#x20;3\. Define  the class, somewhere in the file,  following the template of `ib_walltypes.h`

```cpp
template < typename param, typename cls_t> 
class ibm_user_t;
{

  public:

  compute_surfIB( ...){


  }

};
```
