@testset "Utils" begin
    @testset "Box" begin
        @test_throws ArgumentError box(reshape(Int[], 1, 0))
        @test box(reshape([1], 1, 1)) == [1]
        @test box(reshape([2], 1, 1)) == [1]
        @test box([1 2 4]) == [1,2,4]
        @test box([1 1 2 2; 1 1 2 2]) == [1,1,4,4]
        @test box([1 1 2 2; 1 2 1 2]) == [1,2,3,4]
        @test box([1 2 1 2; 1 1 2 2]) == [1,3,2,4]
        @test box([1 1 1 2 2 2; 1 2 3 1 2 3]) == [1,2,3,4,5,6]
        @test box([1 2 3 1 2 3; 1 1 1 2 2 2]) == [1,3,5,2,4,6]
    end
end
