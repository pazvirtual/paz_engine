const int numEdgeSteps = 10;
const float edgeSteps[] = {1., 1.5, 2., 2., 2., 2., 2., 2., 2., 4.};
const float edgeGuess = 8.;
uniform sampler2D img;
uniform sampler2D lum;
in vec2 uv;
layout(location = 0) out vec4 color;
struct LumData
{
    float m;
    float n;
    float e;
    float s;
    float w;
    float ne;
    float nw;
    float se;
    float sw;
    float minLum;
    float maxLum;
    float contrast;
};
LumData lum_neighborhood(in sampler2D lum, in vec2 texOffset, in vec2 uv)
{
    LumData l;
    l.m = texture(lum, uv).r;
    l.n = texture(lum, uv + texOffset*vec2(0, 1)).r;
    l.e = texture(lum, uv + texOffset*vec2(1, 0)).r;
    l.s = texture(lum, uv + texOffset*vec2(0, -1)).r;
    l.w = texture(lum, uv + texOffset*vec2(-1, 0)).r;
    l.ne = texture(lum, uv + texOffset*vec2(1, 1)).r;
    l.nw = texture(lum, uv + texOffset*vec2(-1, 1)).r;
    l.se = texture(lum, uv + texOffset*vec2(1, -1)).r;
    l.sw = texture(lum, uv + texOffset*vec2(-1, -1)).r;
    l.minLum = min(min(min(min(l.n, l.e), l.s), l.w), l.m);
    l.maxLum = max(max(max(max(l.n, l.e), l.s), l.w), l.m);
    l.contrast = l.maxLum - l.minLum;
    return l;
}
float blend_fac(in LumData l)
{
    float fac = (2.*(l.n + l.e + l.s + l.w) + l.ne + l.nw + l.se + l.sw)/12.;
    fac = abs(fac - l.m);
    fac = clamp(fac/l.contrast, 0., 1.);
    fac = smoothstep(0., 1., fac);
    return fac*fac;
}
struct EdgeData
{
    bool isHorizontal;
    float pixelStep;
    float oppositeLum;
    float grad;
};
EdgeData determine_edge(in sampler2D lum, in vec2 texOffset, in LumData l)
{
    EdgeData e;

    float horizontal = 2.*abs(l.n + l.s - 2.*l.m) + abs(l.ne + l.se - 2.*l.e) +
        abs(l.nw + l.sw - 2.*l.w);
    float vertical = 2.*abs(l.e + l.w - 2.*l.m) + abs(l.ne + l.nw - 2.*l.n) +
        abs(l.se + l.sw - 2.*l.s);
    e.isHorizontal = horizontal >= vertical;

    float pLum = mix(l.e, l.n, float(e.isHorizontal));
    float nLum = mix(l.w, l.s, float(e.isHorizontal));
    float pGrad = abs(pLum - l.m);
    float nGrad = abs(nLum - l.m);
    e.pixelStep = mix(texOffset.x, texOffset.y, float(e.isHorizontal));

    e.pixelStep *= sign(pGrad - nGrad);
    e.oppositeLum = mix(pLum, nLum, float(pGrad < nGrad));
    e.grad = mix(pGrad, nGrad, float(pGrad < nGrad));

    return e;
}
float edge_blend_fac(in sampler2D lum, in vec2 texOffset, in LumData l, in
    EdgeData e, in vec2 uv)
{
    vec2 uvEdge = uv + mix(vec2(0.5*e.pixelStep, 0), vec2(0, 0.5*e.pixelStep),
        float(e.isHorizontal));
    vec2 edgeStep = mix(vec2(0, texOffset.y), vec2(texOffset.x, 0), float(e.
        isHorizontal));

    float edgeLum = 0.5*(l.m + e.oppositeLum);
    float gradThresh = 0.25*e.grad;

    vec2 puv = uvEdge + edgeStep;
    float pLumDelta = texture(lum, puv).r - edgeLum;
    bool pAtEnd = abs(pLumDelta) >= gradThresh;
    for(int i = 0; i + 1 < numEdgeSteps && !pAtEnd; ++i)
    {
        puv += edgeStep*edgeSteps[i];
        pLumDelta = texture(lum, puv).r - edgeLum;
        pAtEnd = abs(pLumDelta) >= gradThresh;
    }
    puv += (1. - float(pAtEnd))*edgeStep*edgeGuess;
    float pDist = mix(puv.y - uv.y, puv.x - uv.x, float(e.isHorizontal));

    vec2 nuv = uvEdge - edgeStep;
    float nLumDelta = texture(lum, nuv).r - edgeLum;
    bool nAtEnd = abs(nLumDelta) >= gradThresh;
    for(int i = 0; i + 1 < numEdgeSteps && !nAtEnd; ++i)
    {
        nuv -= edgeStep*edgeSteps[i];
        nLumDelta = texture(lum, nuv).r - edgeLum;
        nAtEnd = abs(nLumDelta) >= gradThresh;
    }
    nuv -= (1. - float(nAtEnd))*edgeStep*edgeGuess;
    float nDist = mix(uv.y - nuv.y, uv.x - nuv.x, float(e.isHorizontal));

    float shortestDist = mix(nDist, pDist, float(pDist < nDist));
    bool deltaSign = (pDist < nDist && pLumDelta >= 0.) || (pDist >= nDist &&
        nLumDelta >= 0.);

    return float(deltaSign == (l.m - edgeLum >= 0.))*numEdgeSteps*shortestDist;
}
void main()
{
    ivec2 texSize = textureSize(lum, 0);
    vec2 texOffset = vec2(1./texSize.x, 1./texSize.y);
    LumData l = lum_neighborhood(lum, texOffset, uv);
    float pixelFac = blend_fac(l);
    EdgeData e = determine_edge(lum, texOffset, l);
    float edgeFac = edge_blend_fac(lum, texOffset, l, e, uv);
    float fac = max(pixelFac, edgeFac);
    vec2 deltaUv = mix(vec2(e.pixelStep*fac, 0.), vec2(0., e.pixelStep*fac),
        float(e.isHorizontal));
    color = texture(img, uv + deltaUv);
}
