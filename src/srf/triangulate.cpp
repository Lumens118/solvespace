//-----------------------------------------------------------------------------
// Triangulate a surface. If the surface is curved, then we first superimpose
// a grid of quads, with spacing to achieve our chord tolerance. We then
// proceed by ear-clipping; the resulting mesh should be watertight and not
// awful numerically, but has no special properties (Delaunay, etc.).
//
// Copyright 2008-2013 Jonathan Westhues.
//-----------------------------------------------------------------------------
#include "solvespace.h"
#include "gmsh.h"
#include <vector>
#include <list>
#include <fstream>
#include <string>
#include <map>
namespace SolveSpace {

// void ExportContours(const SPolygon& polygon, int loopCount) {
//     FILE* outFile = fopen("d:/work/01Solvespace/SolveSpace02/data1.txt", "a");
//     if(outFile) {
//         fprintf(outFile, "\n========== 循环 %d ==========\n", loopCount);
//         fprintf(outFile, "循环开始时的轮廓数量: %d\n\n", polygon.l.n);

//         int contourIndex = 0;
//         const SContour *c = polygon.l.First();
//         while(c) {
//             fprintf(outFile, "轮廓 %d:\n", contourIndex);
//             fprintf(outFile, "  点数量: %d\n", c->l.n);
//             fprintf(outFile, "  tag: %d\n", c->tag);
//             fprintf(outFile, "  timesEnclosed: %d\n", c->timesEnclosed);

//             Vector currentNormal = c->ComputeNormal();
//             fprintf(outFile, "  法向量: %.6f %.6f %.6f\n",
//                     currentNormal.x, currentNormal.y, currentNormal.z);

//             for(int p = 0; p < c->l.n; p++) {
//                 const Vector& point = c->l[p].p;
//                 fprintf(outFile, "    Point %d: %.6f %.6f %.6f\n",
//                         p, point.x, point.y, point.z);
//             }
//             fprintf(outFile, "\n");

//             c = polygon.l.NextAfter(c);
//             contourIndex++;
//         }
//         fclose(outFile);
//     }
// }

void SPolygon::UvTriangulateInto(SMesh *m, SSurface *srf) {
/*
    思路如下：
    添加三角形到solvespace的方法为AddTriangle()，包含一个STriangle三角元和三个表示三角形xyz顶点的Vector对象。
    从gmsh网格中能获取到的数据结构为点的xyz坐标和每个点对应的索引，还有一个三角形所包括的三个点的索引。
    需要导入进gmsh网格的数据为点集————一个包含点的xyz坐标的数组，和一个或者多个的边界轮廓————边界为一个包含点的索引的数组。
    从solvespace中能够获取到的数据为一个或者多个轮廓————每个轮廓为一个包含点的数组，每个点包括xyz坐标。
*/

//0.初始化    
    // 检查是否要进行处理
    if(l.n <= 0) return;
    // 使用默认的法向量为z轴方向
    normal = Vector::From(0, 0, 1);



//1.获取所有轮廓allContours
//1.1用一个数组保存所有轮廓，按照外层程度排序（越外层索引越小）
    std::vector<SContour*> allContours;

    //计算timeEnclosed————它的意思是轮廓被嵌套的次数，外层为0，内层为1，2，3......
    FixContourDirections();
    //遍历所有轮廓，收集到数组中
    SContour *c = l.First();
    while(c) {
        // 检查轮廓是否有足够的点
        if(c->l.n >= 3) {
            allContours.push_back(c);
        }
        c = l.NextAfter(c);
    }

//1.2获取所有点的坐标和所有轮廓的边
    std::vector<double> coordinates3d;//x0,y0,z0;x1,y1,z1......所有点的坐标
    std::vector<int> outerEdges3d;//外轮廓的边（逆时针）
    std::vector<std::vector<int>> innerEdges3d;//内轮廓的边（顺时针），每个元素是一个内轮廓的边数组edges3d（其实内外轮廓不需要方向，之前误以为需要）
    std::vector<SContour*> innerContour3d;//存储内轮廓的数组
    std::vector<SContour*> outerContour3d;//存储外轮廓的数组

//1.2.1 区分内轮廓和外轮廓innerContour3d和outerContour3d
    for(size_t i = 0; i < allContours.size(); i++) {
        SContour* contour = allContours[i];
        if(contour->timesEnclosed == 0) {
            outerContour3d.push_back(contour);
        } else {
            innerContour3d.push_back(contour);
        }
    }

//1.2.2 把轮廓上所有点的坐标全部传进coordinates3d
    // 临时的比较工具，用于比较Vector对象
    struct VectorCompare {
        bool operator()(const Vector& a, const Vector& b) const {
            // 按x, y, z坐标依次比较
            if(a.x != b.x) return a.x < b.x;
            if(a.y != b.y) return a.y < b.y;
            return a.z < b.z;
        }
    };
    std::map<Vector, std::size_t, VectorCompare> pointIndexMap;
    std::size_t currentIndex = 1;

    // 先处理外轮廓，将所有点的坐标添加到coordinates3d
    for(size_t i = 0; i < outerContour3d.size(); i++) {
        SContour* contour = outerContour3d[i];
        // 保存点信息
        for(int j = 0; j < contour->l.n; j++) {
            const Vector& point = contour->l[j].p;
            // 检查点是否已存在
            if(pointIndexMap.find(point) == pointIndexMap.end()) {
                // 添加新点到coordinates3d
                coordinates3d.push_back(point.x);
                coordinates3d.push_back(point.y);
                coordinates3d.push_back(point.z);
                // 记录点索引（从1开始）
                pointIndexMap[point] = currentIndex++;
            }
        }

        // 构建外轮廓的边，确保轮廓闭合且逆时针
        for(int j = 0; j < contour->l.n; j++) {
            const Vector& currentPoint = contour->l[j].p;
            // 下一个点，最后一个点连接回第一个点
            int nextIndex = (j + 1) % contour->l.n;
            const Vector& nextPoint = contour->l[nextIndex].p;

            // 确保两个点都在索引映射中
            if(pointIndexMap.find(currentPoint) != pointIndexMap.end() && 
               pointIndexMap.find(nextPoint) != pointIndexMap.end()) {
                // 添加外轮廓边（点索引对）到outerEdges3d（保持原逆时针方向）
                outerEdges3d.push_back(static_cast<int>(pointIndexMap[currentPoint]));
                outerEdges3d.push_back(static_cast<int>(pointIndexMap[nextPoint]));
            }
        }
    }

    // 处理内轮廓，先添加点，再构建边
    for(size_t innerIdx = 0; innerIdx < innerContour3d.size(); innerIdx++) {
        SContour* contour = innerContour3d[innerIdx];

        // 先添加内轮廓的点到coordinates3d
        for(int j = 0; j < contour->l.n; j++) {
            const Vector& point = contour->l[j].p;
            // 检查点是否已存在
            if(pointIndexMap.find(point) == pointIndexMap.end()) {
                // 添加新点到coordinates3d
                coordinates3d.push_back(point.x);
                coordinates3d.push_back(point.y);
                coordinates3d.push_back(point.z);
                // 记录点索引（从1开始）
                pointIndexMap[point] = currentIndex++;
            }
        }
    }

    // 再处理内轮廓，构建innerEdges3d（确保顺时针方向）
    for(size_t innerIdx = 0; innerIdx < innerContour3d.size(); innerIdx++) {
        SContour* contour = innerContour3d[innerIdx];
        std::vector<int> currentInnerEdges;

        // 构建内轮廓的边，先按原顺序添加
        for(int j = 0; j < contour->l.n; j++) {
            const Vector& currentPoint = contour->l[j].p;
            // 下一个点，最后一个点连接回第一个点
            int nextIndex = (j + 1) % contour->l.n;
            const Vector& nextPoint = contour->l[nextIndex].p;

            // 确保两个点都在索引映射中
            if(pointIndexMap.find(currentPoint) != pointIndexMap.end() && 
               pointIndexMap.find(nextPoint) != pointIndexMap.end()) {
                currentInnerEdges.push_back(static_cast<int>(pointIndexMap[currentPoint]));
                currentInnerEdges.push_back(static_cast<int>(pointIndexMap[nextPoint]));
            }
        }

        // 反转内轮廓的边顺序，确保顺时针方向
        std::vector<int> reversedEdges;
        for(int k = static_cast<int>(currentInnerEdges.size()) - 2; k >= 0; k -= 2) {
            // 反转每条边的起点和终点
            reversedEdges.push_back(currentInnerEdges[k + 1]);
            reversedEdges.push_back(currentInnerEdges[k]);
        }

        // 将当前内轮廓的边添加到innerEdges3d
        innerEdges3d.push_back(reversedEdges);
    }



    //2. 用gmsh进行处理
    //2.1 初始化gmsh
    if(!gmsh::isInitialized()){
        gmsh::initialize();
    }

    //2.2 创建模型
    gmsh::model::add("Triangulation");

    //设置特征长度characteristic length
    double lc = 0.1;

    //2.2.1 创建点points
    std::vector<std::vector<double>> points;
    //把coordinates3d里的内容转移到points里
    for(size_t j = 0; j < coordinates3d.size(); j += 3) {
        points.push_back({coordinates3d[j], coordinates3d[j+1], coordinates3d[j+2]});
    }
    //添加点到gmsh
    for(size_t j = 0; j < points.size(); j++) {
        gmsh::model::geo::addPoint(points[j][0], points[j][1], points[j][2], lc, j+1);
    }
    std::map<std::pair<int, int>, int> curveTagMap;  // 存储点对到曲线标签的映射
    int currentCurveTag = 1;

    //2.2.2 创建外轮廓的曲线outerCurves
    std::vector<int> outerCurves;
    for(size_t j = 0; j < outerEdges3d.size(); j += 2) {
        int startPoint = outerEdges3d[j];
        int endPoint = outerEdges3d[j+1];

        // 创建曲线（边）
        int curveTag = gmsh::model::geo::addLine(startPoint, endPoint, currentCurveTag++);
        outerCurves.push_back(curveTag);
        curveTagMap[{startPoint, endPoint}] = curveTag;
    }
    // 创建外曲线环
    gmsh::model::geo::addCurveLoop(outerCurves, 1); 

    //2.2.3 创建内轮廓的曲线innerCurves
    std::vector<int> innerCurves;
    for(size_t j = 0; j < innerEdges3d.size(); j++) {
        std::vector<int> currentCurves;
        for(size_t k = 0; k < innerEdges3d[j].size(); k += 2) {
            int startPoint = innerEdges3d[j][k];
            int endPoint = innerEdges3d[j][k+1];

            // 创建曲线（边）
            int curveTag = gmsh::model::geo::addLine(startPoint, endPoint, currentCurveTag++);
            currentCurves.push_back(curveTag);
        }

        // 创建内曲线环
        int innerLoopTag = gmsh::model::geo::addCurveLoop(currentCurves, j+2);
        innerCurves.push_back(innerLoopTag);
    }

    //2.2.4 创建平面
    std::vector<int> allCurveLoops;
    // 外边界曲线环标签
    allCurveLoops.push_back(1);
    // 添加所有内边界曲线环标签
    for(int loopTag : innerCurves) {
        allCurveLoops.push_back(loopTag);
    }
    // 创建包含内外边界的平面域
    gmsh::model::geo::addPlaneSurface(allCurveLoops, 1);

    //2.2.5 同步几何
    gmsh::model::geo::synchronize();

    //2.3 处理网格
    gmsh::option::setNumber("Mesh.Algorithm", 5);  // 5: Delaunay算法
    gmsh::model::mesh::generate(2);

    //2.4 提取三角剖分网格的数据
    std::vector<std::size_t> nodeTags;
    std::vector<double> nodeCoords;
    std::vector<double> parametricCoords;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoords, parametricCoords, -1, -1, false, false);

    // 创建节点标签到索引的映射
    std::map<std::size_t, int> nodeMap;
    for(int i = 0; i < nodeTags.size(); i++) {
        nodeMap[nodeTags[i]] = i;
    }

    // 获取所有元素
    std::vector<int> elementTypes;
    std::vector<std::vector<std::size_t>> elementTags;
    std::vector<std::vector<std::size_t>> elementNodeTags;
    gmsh::model::mesh::getElements(elementTypes, elementTags, elementNodeTags, -1, -1);

    //2.5 三角形元素（类型为2）
    for(size_t i = 0; i < elementTypes.size(); i++) {
        if(elementTypes[i] == 2) { // 三角形元素类型为2
            const std::vector<std::size_t>& triangleNodeTags = elementNodeTags[i];

            // 每个三角形有3个节点
            for(int j = 0; j < triangleNodeTags.size(); j += 3) {
                // 获取三角形的三个节点索引
                std::size_t tag1 = triangleNodeTags[j];
                std::size_t tag2 = triangleNodeTags[j+1];
                std::size_t tag3 = triangleNodeTags[j+2];

                // 从节点标签映射中获取对应的坐标索引
                int idx1 = nodeMap[tag1];
                int idx2 = nodeMap[tag2];
                int idx3 = nodeMap[tag3];

                // 创建三个顶点
                Vector v1 = Vector::From(nodeCoords[3*idx1], nodeCoords[3*idx1+1], nodeCoords[3*idx1+2]);
                Vector v2 = Vector::From(nodeCoords[3*idx2], nodeCoords[3*idx2+1], nodeCoords[3*idx2+2]);
                Vector v3 = Vector::From(nodeCoords[3*idx3], nodeCoords[3*idx3+1], nodeCoords[3*idx3+2]);

                // 添加三角形到网格，使用表面的颜色
                STriMeta meta = { srf->face, srf->color };
                m->AddTriangle(meta, v1, v2, v3);
            }
        }
    }

    //2.6 清理Gmsh模型
    gmsh::model::remove();
    gmsh::finalize();
}


    


bool SContour::BridgeToContour(SContour *sc,
                               SEdgeList *avoidEdges, List<Vector> *avoidPts)
{
    int i, j;
    bool withbridge = true;

    // Start looking for a bridge on our new hole near its leftmost (min x)
    // point.
    int sco = 0;
    for(i = 0; i < (sc->l.n - 1); i++) {
        if((sc->l[i].p).EqualsExactly(sc->xminPt)) {
            sco = i;
        }
    }

    // And start looking on our merged contour at whichever point is nearest
    // to the leftmost point of the new segment.
    int thiso = 0;
    double dmin = 1e10;
    for(i = 0; i < l.n-1; i++) {
        Vector p = l[i].p;
        double d = (p.Minus(sc->xminPt)).MagSquared();
        if(d < dmin) {
            dmin = d;
            thiso = i;
        }
    }

    int thisp, scp;

    Vector a, b, *f;

    // First check if the contours share a point; in that case we should
    // merge them there, without a bridge.
    for(i = 0; i < l.n; i++) {
        thisp = WRAP(i+thiso, l.n-1);
        a = l[thisp].p;

        for(f = avoidPts->First(); f; f = avoidPts->NextAfter(f)) {
            if(f->Equals(a)) break;
        }
        if(f) continue;

        for(j = 0; j < (sc->l.n - 1); j++) {
            scp = WRAP(j+sco, (sc->l.n - 1));
            b = sc->l[scp].p;

            if(a.Equals(b)) {
                withbridge = false;
                goto haveEdge;
            }
        }
    }

    // If that fails, look for a bridge that does not intersect any edges.
    for(i = 0; i < l.n; i++) {
        thisp = WRAP(i+thiso, l.n);
        a = l[thisp].p;

        for(f = avoidPts->First(); f; f = avoidPts->NextAfter(f)) {
            if(f->Equals(a)) break;
        }
        if(f) continue;

        for(j = 0; j < (sc->l.n - 1); j++) {
            scp = WRAP(j+sco, (sc->l.n - 1));
            b = sc->l[scp].p;

            for(f = avoidPts->First(); f; f = avoidPts->NextAfter(f)) {
                if(f->Equals(b)) break;
            }
            if(f) continue;

            if(avoidEdges->AnyEdgeCrossings(a, b) > 0) {
                // doesn't work, bridge crosses an existing edge
            } else {
                goto haveEdge;
            }
        }
    }

    // Tried all the possibilities, didn't find an edge
    return false;

haveEdge:
    SContour merged = {};
    for(i = 0; i < l.n; i++) {
        if(withbridge || (i != thisp)) {
            merged.AddPoint(l[i].p);
        }
        if(i == thisp) {
            // less than or equal; need to duplicate the join point
            for(j = 0; j <= (sc->l.n - 1); j++) {
                int jp = WRAP(j + scp, (sc->l.n - 1));
                merged.AddPoint((sc->l[jp]).p);
            }
            // and likewise duplicate join point for the outer curve
            if(withbridge) {
                merged.AddPoint(l[i].p);
            }
        }
    }

    // and future bridges mustn't cross our bridge, and it's tricky to get
    // things right if two bridges come from the same point
    if(withbridge) {
        avoidEdges->AddEdge(a, b);
        avoidPts->Add(&a);
    }
    avoidPts->Add(&b);

    l.Clear();
    l = merged.l;
    return true;
}

bool SContour::IsEmptyTriangle(int ap, int bp, int cp, double scaledEPS) const {

    STriangle tr = {};
    tr.a = l[ap].p;
    tr.b = l[bp].p;
    tr.c = l[cp].p;

    // Accelerate with an axis-aligned bounding box test
    Vector maxv = tr.a, minv = tr.a;
    (tr.b).MakeMaxMin(&maxv, &minv);
    (tr.c).MakeMaxMin(&maxv, &minv);

    Vector n = Vector::From(0, 0, -1);

    int i;
    for(i = 0; i < l.n; i++) {
        if(i == ap || i == bp || i == cp) continue;

        Vector p = l[i].p;
        if(p.OutsideAndNotOn(maxv, minv)) continue;

        // A point on the edge of the triangle is considered to be inside,
        // and therefore makes it a non-ear; but a point on the vertex is
        // "outside", since that's necessary to make bridges work.
        if(p.EqualsExactly(tr.a)) continue;
        if(p.EqualsExactly(tr.b)) continue;
        if(p.EqualsExactly(tr.c)) continue;

        if(tr.ContainsPointProjd(n, p)) {
            return false;
        }
    }
    return true;
}

// Test if ray b->d passes through triangle a,b,c
static bool RayIsInside(Vector a, Vector c, Vector b, Vector d) {
    // coincident edges are not considered to intersect the triangle
    if (d.Equals(a)) return false;
    if (d.Equals(c)) return false;
    // if d and c are on opposite sides of ba, we are ok
    // likewise if d and a are on opposite sides of bc
    Vector ba = a.Minus(b);
    Vector bc = c.Minus(b);
    Vector bd = d.Minus(b);

    // perpendicular to (x,y) is (x,-y) so dot that with the two points. If they
    // have opposite signs their product will be negative. If bd and bc are on
    // opposite sides of ba the ray does not intersect. Likewise for bd,ba and bc.
    if ( (bd.x*(ba.y) + (bd.y * (-ba.x))) * ( bc.x*(ba.y) + (bc.y * (-ba.x))) < LENGTH_EPS)
        return false;
    if ( (bd.x*(bc.y) + (bd.y * (-bc.x))) * ( ba.x*(bc.y) + (ba.y * (-bc.x))) < LENGTH_EPS)
        return false;

    return true;
}

// bool SContour::IsEar(int bp, double scaledEps) const {
//     int ap = WRAP(bp-1, l.n),
//         cp = WRAP(bp+1, l.n);

//     STriangle tr = {};
//     tr.a = l[ap].p;
//     tr.b = l[bp].p;
//     tr.c = l[cp].p;

//     if((tr.a).Equals(tr.c)) {
//         // This is two coincident and anti-parallel edges. Zero-area, so
//         // won't generate a real triangle, but we certainly can clip it.
//         return true;
//     }

//     Vector n = Vector::From(0, 0, -1);
//     if((tr.Normal()).Dot(n) < scaledEps) {
//         // This vertex is reflex, or between two collinear edges; either way,
//         // it's not an ear.
//         return false;
//     }

//     // Accelerate with an axis-aligned bounding box test
//     Vector maxv = tr.a, minv = tr.a;
//     (tr.b).MakeMaxMin(&maxv, &minv);
//     (tr.c).MakeMaxMin(&maxv, &minv);

//     int i;
//     for(i = 0; i < l.n; i++) {
//         if(i == ap || i == bp || i == cp) continue;

//         Vector p = l[i].p;
//         if(p.OutsideAndNotOn(maxv, minv)) continue;

//         // A point on the edge of the triangle is considered to be inside,
//         // and therefore makes it a non-ear; but a point on the vertex is
//         // "outside", since that's necessary to make bridges work.
//         if(p.EqualsExactly(tr.a)) continue;
//         if(p.EqualsExactly(tr.c)) continue;
//         // points coincident with bp have to be allowed for bridges but edges
//         // from that other point must not cross through our triangle.
//         if(p.EqualsExactly(tr.b)) {
//             int j = WRAP(i-1, l.n);
//             int k = WRAP(i+1, l.n);
//             Vector jp = l[j].p;
//             Vector kp = l[k].p;

//             // two consecutive bridges (A,B,C) and later (C,B,A) are not an ear
//             if (jp.Equals(tr.c) && kp.Equals(tr.a)) return false;
//             // check both edges from the point in question
//             if (!RayIsInside(tr.a, tr.c, p,jp) && !RayIsInside(tr.a, tr.c, p,kp))
//                 continue;
//         }

//         if(tr.ContainsPointProjd(n, p)) {
//             return false;
//         }
//     }
//     return true;
// }

void SContour::ClipEarInto(SMesh *m, int bp, double scaledEps) {
    int ap = WRAP(bp-1, l.n),
        cp = WRAP(bp+1, l.n);

    STriangle tr = {};
    tr.a = l[ap].p;
    tr.b = l[bp].p;
    tr.c = l[cp].p;
    if(tr.Normal().MagSquared() < scaledEps*scaledEps) {
        // A vertex with more than two edges will cause us to generate
        // zero-area triangles, which must be culled.
    } else {
        m->AddTriangle(&tr);
    }

    // By deleting the point at bp, we may change the ear-ness of the points
    // on either side.
    l[ap].ear = EarType::UNKNOWN;
    l[cp].ear = EarType::UNKNOWN;

    l.ClearTags();
    l[bp].tag = 1;
    l.RemoveTagged();
}

// void SContour::UvTriangulateInto(SMesh *m, SSurface *srf) {
//     Vector tu, tv;
//     srf->TangentsAt(0.5, 0.5, &tu, &tv);
//     double s = sqrt(tu.MagSquared() + tv.MagSquared());
//     // We would like to apply our tolerances in xyz; but that would be a lot
//     // of work, so at least scale the epsilon semi-reasonably. That's
//     // perfect for square planes, less perfect for anything else.
//     double scaledEps = LENGTH_EPS / s;

//     int i;
//     // Clean the original contour by removing any zero-length edges.
//     // initialize eartypes to unknown while we're going over them.
//     l.ClearTags();
//     l[0].ear = EarType::UNKNOWN;
//     for(i = 1; i < l.n; i++) {
//        l[i].ear = EarType::UNKNOWN;
//        if((l[i].p).Equals(l[i-1].p)) {
//             l[i].tag = 1;
//         }
//     }
//     if( (l[0].p).Equals(l[l.n-1].p) ) {
//         l[l.n-1].tag = 1;
//     }
//     l.RemoveTagged();

//     // Handle simple triangle fans all at once. This pass is optional.
//     if(srf->degm == 1 && srf->degn == 1) {
//         l.ClearTags();
//         int j=0;
//         int pstart = 0;
//         double elen = -1.0;
//         double oldspan = 0.0;
//         for(i = 1; i < l.n; i++) {
//             Vector ab = l[i].p.Minus(l[i-1].p);
//             // first time just measure the segment
//             if (elen < 0.0) {
//                 elen = ab.Dot(ab);
//                 oldspan = elen;
//                 j = 1;
//                 continue;
//             }
//             // check for consecutive segments of similar size which are also
//             // ears and where the group forms a convex ear
//             bool end = false;
//             double ratio = ab.Dot(ab) / elen;
//             if ((ratio < 0.25) || (ratio > 4.0)) end = true;

//             double slen = l[pstart].p.Minus(l[i].p).MagSquared();
//             if (slen < oldspan) end = true;

//             if (!IsEar(i-1, scaledEps) ) end = true;
// //            if ((j>0) && !IsEar(pstart, i-1, i, scaledEps)) end = true;
//             if ((j>0) && !IsEmptyTriangle(pstart, i-1, i, scaledEps)) end = true;
//             // the new segment is valid so add to the fan
//             if (!end) {
//                 j++;
//                 oldspan = slen;
//             }
//             // we need to stop at the end of polygon but may still
//             if (i == l.n-1) {
//                 end = true;
//             }
//             if (end) {  // triangulate the fan and tag the vertices
//                 if (j > 3) {
//                     Vector center = l[pstart+1].p.Plus(l[pstart+j-1].p).ScaledBy(0.5);
//                     for (int x=0; x<j; x++) {
//                         STriangle tr = {};
//                         tr.a = center;
//                         tr.b = l[pstart+x].p;
//                         tr.c = l[pstart+x+1].p;
//                         m->AddTriangle(&tr);
//                     }
//                     for (int x=1; x<j; x++) {
//                         l[pstart+x].tag = 1;
//                     }
//                     STriangle tr = {};
//                     tr.a = center;
//                     tr.b = l[pstart+j].p;
//                     tr.c = l[pstart].p;
//                     m->AddTriangle(&tr);    
//                 }
//                 pstart = i-1;
//                 elen = ab.Dot(ab);
//                 oldspan = elen;
//                 j = 1;
//             }
//         }
//         l.RemoveTagged();
//     }  // end optional fan creation pass

//     bool toggle = false;
//     while(l.n > 3) {
//         int bestEar = -1;
//         double bestChordTol = VERY_POSITIVE;
//         // Alternate the starting position so we generate strip-like
//         // triangulations instead of fan-like
//         toggle = !toggle;
//         int offset = toggle ? -1 : 0;
//         for(i = 0; i < l.n; i++) {
//             int ear = WRAP(i+offset, l.n);
//             if(l[ear].ear == EarType::UNKNOWN) {
//                 (l[ear]).ear = IsEar(ear, scaledEps) ? EarType::EAR : EarType::NOT_EAR;
//             }
//             if(l[ear].ear == EarType::EAR) {
//                 if(srf->degm == 1 && srf->degn == 1) {
//                     // This is a plane; any ear is a good ear.
//                     bestEar = ear;
//                     break;
//                 }
//                 // If we are triangulating a curved surface, then try to
//                 // clip ears that have a small chord tolerance from the
//                 // surface.
//                 Vector prev = l[WRAP((i+offset-1), l.n)].p,
//                        next = l[WRAP((i+offset+1), l.n)].p;
//                 double tol = srf->ChordToleranceForEdge(prev, next);
//                 if(tol < bestChordTol - scaledEps) {
//                     bestEar = ear;
//                     bestChordTol = tol;
//                 }
//                 if(bestChordTol < 0.1*SS.ChordTolMm()) {
//                     break;
//                 }
//             }
//         }
//         if(bestEar < 0) {
//             dbp("couldn't find an ear! fail");
//             return;
//         }
//         ClipEarInto(m, bestEar, scaledEps);
//     }

//     ClipEarInto(m, 0, scaledEps); // add the last triangle
// }

double SSurface::ChordToleranceForEdge(Vector a, Vector b) const {
    Vector as = PointAt(a.x, a.y), bs = PointAt(b.x, b.y);

    double worst = VERY_NEGATIVE;
    int i;
    for(i = 1; i <= 3; i++) {
        Vector p  = a. Plus((b. Minus(a )).ScaledBy(i/4.0)),
               ps = as.Plus((bs.Minus(as)).ScaledBy(i/4.0));

        Vector pps = PointAt(p.x, p.y);
        worst = max(worst, (pps.Minus(ps)).MagSquared());
    }
    return sqrt(worst);
}

Vector SSurface::PointAtMaybeSwapped(double u, double v, bool swapped) const {
    if(swapped) {
        return PointAt(v, u);
    } else {
        return PointAt(u, v);
    }
}

Vector SSurface::NormalAtMaybeSwapped(double u, double v, bool swapped) const {
    Vector du, dv;
    if(swapped) {
        TangentsAt(v, u, &dv, &du);
    } else {
        TangentsAt(u, v, &du, &dv);
    }
    return du.Cross(dv).WithMagnitude(1.0);
}

void SSurface::MakeTriangulationGridInto(List<double> *l, double vs, double vf,
                                         bool swapped, int depth) const
{
    double worst = 0;

    // Try piecewise linearizing four curves, at u = 0, 1/3, 2/3, 1; choose
    // the worst chord tolerance of any of those.
    double worst_twist = 1.0;
    int i;
    for(i = 0; i <= 3; i++) {
        double u = i/3.0;

        // This chord test should be identical to the one in SBezier::MakePwl
        // to make the piecewise linear edges line up with the grid more or
        // less.
        Vector ps = PointAtMaybeSwapped(u, vs, swapped),
               pf = PointAtMaybeSwapped(u, vf, swapped);

        double vm1 = (2*vs + vf) / 3,
               vm2 = (vs + 2*vf) / 3;

        Vector pm1 = PointAtMaybeSwapped(u, vm1, swapped),
               pm2 = PointAtMaybeSwapped(u, vm2, swapped);

        // 0.999 is about 2.5 degrees of twist over the middle 1/3 V-span.
        // we don't check at the ends because the derivative may not be valid there.
        double twist = 1.0;
        if (degm == 1) twist = NormalAtMaybeSwapped(u, vm1, swapped).Dot(
                               NormalAtMaybeSwapped(u, vm2, swapped) );
        if (twist < worst_twist) worst_twist = twist;

        worst = max(worst, pm1.DistanceToLine(ps, pf.Minus(ps)));
        worst = max(worst, pm2.DistanceToLine(ps, pf.Minus(ps)));
    }

    double step = 1.0/SS.GetMaxSegments();
    if( ((vf - vs) < step || worst < SS.ChordTolMm())
        && ((worst_twist > 0.999) || (depth > 3)) ) {
        l->Add(&vf);
    } else {
        MakeTriangulationGridInto(l, vs, (vs+vf)/2, swapped, depth+1);
        MakeTriangulationGridInto(l, (vs+vf)/2, vf, swapped, depth+1);
    }
}

void SPolygon::UvGridTriangulateInto(SMesh *mesh, SSurface *srf) {
    SEdgeList orig = {};
    MakeEdgesInto(&orig);

    SEdgeList holes = {};

    normal = Vector::From(0, 0, 1);
    FixContourDirections();

    // Build a rectangular grid, with horizontal and vertical lines in the
    // uv plane. The spacing of these lines is adaptive, so calculate that.
    List<double> li, lj;
    li = {};
    lj = {};
    double v[5] = {0.0, 0.25, 0.5, 0.75, 1.0};
    li.Add(&v[0]);
    srf->MakeTriangulationGridInto(&li, 0, 1, /*swapped=*/true, 0);
    lj.Add(&v[0]);
    srf->MakeTriangulationGridInto(&lj, 0, 1, /*swapped=*/false, 0);

    // force 2nd order grid to have at least 4 segments in each direction
    if ((li.n < 5) && (srf->degm>1)) { // 4 segments minimum
        li.Clear();
        li.Add(&v[0]);li.Add(&v[1]);li.Add(&v[2]);li.Add(&v[3]);li.Add(&v[4]);
    }
    if ((lj.n < 5) && (srf->degn>1)) { // 4 segments minimum
        lj.Clear();
        lj.Add(&v[0]);lj.Add(&v[1]);lj.Add(&v[2]);lj.Add(&v[3]);lj.Add(&v[4]);
    }

    if ((li.n > 3) && (lj.n > 3)) {
        // Now iterate over each quad in the grid. If it's outside the polygon,
        // or if it intersects the polygon, then we discard it. Otherwise we
        // generate two triangles in the mesh, and cut it out of our polygon.
        // Quads around the perimeter would be rejected by AnyEdgeCrossings.
        std::vector<bool> bottom(lj.n, false); // did we use this quad?
        Vector tu = {0,0,0}, tv = {0,0,0};
        int i, j;
        for(i = 1; i < (li.n-1); i++) {
            bool prev_flag = false;
            for(j = 1; j < (lj.n-1); j++) {
                bool this_flag = true;
                double us = li[i], uf = li[i+1],
                       vs = lj[j], vf = lj[j+1];

                Vector a = Vector::From(us, vs, 0),
                       b = Vector::From(us, vf, 0),
                       c = Vector::From(uf, vf, 0),
                       d = Vector::From(uf, vs, 0);

                //  |   d-----c
                //  |   |     |
                //  |   |     |
                //  |   a-----b
                //  |
                //  +-------------> j/v axis

                if( (i==(li.n-2)) || (j==(lj.n-2)) ||
                   orig.AnyEdgeCrossings(a, b, NULL) ||
                   orig.AnyEdgeCrossings(b, c, NULL) ||
                   orig.AnyEdgeCrossings(c, d, NULL) ||
                   orig.AnyEdgeCrossings(d, a, NULL))
                {
                    this_flag = false;
                }

                // There's no intersections, so it doesn't matter which point
                // we decide to test.
                if(!this->ContainsPoint(a)) {
                    this_flag = false;
                }
                
                if (this_flag) {
                    // Add the quad to our mesh
                    srf->TangentsAt(us,vs, &tu,&tv);
                    if(tu.Dot(tv) < LENGTH_EPS) {
                        /* Split "the other way" if angle>90
                           compare to LENGTH_EPS instead of zero to avoid alternating triangle
                           "orientations" when the tangents are orthogonal (revolve, lathe etc.)
                           this results in a higher quality mesh. */
                        STriangle tr = {};
                        tr.a = a;
                        tr.b = b;
                        tr.c = c;
                        mesh->AddTriangle(&tr);
                        tr.a = a;
                        tr.b = c;
                        tr.c = d;
                        mesh->AddTriangle(&tr);
                    } else{
                        STriangle tr = {};
                        tr.a = a;
                        tr.b = b;
                        tr.c = d;
                        mesh->AddTriangle(&tr);
                        tr.a = b;
                        tr.b = c;
                        tr.c = d;
                        mesh->AddTriangle(&tr);
                    }
                    if (!prev_flag) // add our own left edge
                        holes.AddEdge(d, a);
                    if (!bottom[j]) // add our own bottom edge
                        holes.AddEdge(a, b);
                } else {
                    if (prev_flag) // add our left neighbors right edge
                        holes.AddEdge(a, d);
                    if (bottom[j]) // add our bottom neighbors top edge
                        holes.AddEdge(b, a);
                }
                prev_flag = this_flag;
                bottom[j] = this_flag;
            }
        }

        // Because no duplicate edges were created we do not need to cull them.
        SPolygon hp = {};
        holes.AssemblePolygon(&hp, NULL, /*keepDir=*/true);

        SContour *sc;
        for(sc = hp.l.First(); sc; sc = hp.l.NextAfter(sc)) {
            l.Add(sc);
        }
        hp.l.Clear();
    }
    orig.Clear();
    holes.Clear();
    li.Clear();
    lj.Clear();
    UvTriangulateInto(mesh, srf);
}

void SPolygon::TriangulateInto(SMesh *m) const {
    Vector n = normal;
    if(n.Equals(Vector::From(0.0, 0.0, 0.0))) {
       n = ComputeNormal();
    }
    Vector u = n.Normal(0);
    Vector v = n.Normal(1);

    SPolygon p = {};
    this->InverseTransformInto(&p, u, v, n);

    SSurface srf = SSurface::FromPlane(Vector::From(0.0, 0.0, 0.0),
                                       Vector::From(1.0, 0.0, 0.0),
                                       Vector::From(0.0, 1.0, 0.0));
    SMesh pm = {};
    p.UvTriangulateInto(&pm, &srf);
    for(STriangle st : pm.l) {
        st = st.Transform(u, v, n);
        m->AddTriangle(&st);
    }

    p.Clear();
    pm.Clear();
}

} // namespace SolveSpace










