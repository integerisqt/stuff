local err = 10e-9;
local acos = math.acos;
local cos = math.cos;

local function iszero(x)
    return x > -err and x < err;
end

local function cubesqrt(x)
    return x ^ ((1 / 3) * (x > 0 and 1 or -1));
end

local function solvequadric(c0, c1, c2)
    local s0, s1;
    local p, q, D;
    p = c1 / (2 * c0);
    q = c2 / c0;
    D = p * p - q;
    if iszero(D) then
        s0 = -p;
        return 0;
    elseif D < 0 then
        return;
    else
        local sqrt_D = D ^ 0.5;
        s0 = sqrt_D - p;
        s1 = -sqrt_D - p;
        return s0, s1;
    end
end

local function solvequbic(c0, c1, c2, c3)
    local s0, s1, s2;
    local num, sub;
    local A, B, C;
    local sq_A, p, q;
    local cb_p, D;
    A = c1 / c0;
    B = c2 / c0;
    C = c3 / c0;

    sq_A = A * A;
    p = (1 / 3) * (-(1 / 3) * sq_A + B);
    q = 0.5 * ((2 / 27) * A * sq_A - (1 / 3) * A * B + C);

    cb_p = p * p * p;
    D = q * q + cb_p;

    if iszero(D) then
        if iszero(q) then
            s0 = 0;
            num = 1;
        else
            local u = cubesqrt(-q);
            s0 = 2 * u;
            s1 = -u;
            num = 2;
        end
    elseif D < 0 then
        local phi = (1 / 3) * acos(-q / ((-cb_p) ^ 0.5));
        local t = 2 * ((-p) ^ 0.5);

        s0 = t * cos(phi);
        s1 = -t * cos(phi + math.pi / 3);
        s2 = -t * cos(phi - math.pi / 3);
        num = 3;
    else
        local sqrt_D = D ^ 0.5;
        local u = cubesqrt(sqrt_D - q);
        local v = -cubesqrt(sqrt_D + q);
        s0 = u + v;
        num = 1;
    end
    sub = (1 / 3) * A;
    if num > 0 then s0 = s0 - sub; end
    if num > 1 then s1 = s1 - sub; end
    if num > 2 then s2 = s2 - sub; end
    return s0, s1, s2;
end

local function solvequartic(c0, c1, c2, c3 ,c4)
    local s0, s1, s2, s3;

    local coeffs = {};
    local z, u, v, sub;
    local A, B, C, D;
    local sq_A, p, q, r;
    local num;
    A = c1 / c0;
    B = c2 / c0;
    C = c3 / c0;
    D = c4 / c0;

    sq_A = A * A;
    p = -0.375 * sq_A + B;
    q = 0.125 * sq_A * A - 0.5 * A * B + C;
    r = -(3 / 256) * sq_A * sq_A + 0.0625 * sq_A * B - 0.25 * A * C + D;

    if iszero(r) then
        coeffs[3] = q;
        coeffs[2] = p;
        coeffs[1] = 0;
        coeffs[0] = 1;
        local results = {solvequbic(coeffs[0], coeffs[1], coeffs[2], coeffs[3])};
        num = #results;
        s0, s1, s2 = results[1], results[2], results[3];
    else
        coeffs[3] = 0.5 * r * p - 0.125 * q * q;
        coeffs[2] = -r;
        coeffs[1] = -0.5 * p;
        coeffs[0] = 1;
        s0, s1, s2 = solvequbic(coeffs[0], coeffs[1], coeffs[2], coeffs[3]);
        z = s0;
        u = z * z - r;
        v = 2 * z - p;
        if iszero(u) then
            u = 0;
        elseif u > 0 then
            u = u ^ 0.5;
        else
            return;
        end
        if iszero(v) then
            v = 0;
        elseif v > 0 then
            v = v ^ 0.5;
        else
            return;
        end
        coeffs[2] = z - u;
        coeffs[1] = q < 0 and -v or v;
        coeffs[0] = 1;

        do
            local results = {solvequadric(coeffs[0], coeffs[1], coeffs[2])};
            num = #results;
            s0, s1 = results[1], results[2];
        end

        coeffs[2] = z + u;
        coeffs[1] = q < 0 and v or -v;
        coeffs[0] = 1;

        if num == 0 then
            local results = {solvequadric(coeffs[0], coeffs[1], coeffs[2])};
            num = num + #results;
            s0, s1 = results[1], results[2];
        end
        if num == 1 then
            local results = {solvequadric(coeffs[0], coeffs[1], coeffs[2])};
            num = num + #results;
            s1, s2 = results[1], results[2];
        end
        if num == 2 then
            local results = {solvequadric(coeffs[0], coeffs[1], coeffs[2])};
            num = num + #results;
            s2, s3 = results[1], results[2];
        end
    end
    sub = 0.25 * A;
    if num > 0 then s0 = s0 - sub; end
	if num > 1 then s1 = s1 - sub; end
	if num > 2 then s2 = s2 - sub; end
	if num > 3 then s3 = s3 - sub; end

	return s0, s1, s2, s3;
end

return solvequartic;
