#include <iostream>
#include "gtest/gtest.h"
#include <math.h>
#include <vector>
#include <random>
#include <utils/PCA.hpp>


namespace {
	TEST(Directions, SameDirections) {
		OnlinePCA pca;
		
		glm::vec3 vec(1.f, 1.f, 1.f);
		glm::vec3 norm = glm::normalize(vec);
		pca.addVector(vec);
		glm::vec3 direction = pca.getPrincipalDirection();
		EXPECT_EQ(direction, norm);	// first direction should be equal to normalized vector with no ambiguity

		vec = glm::vec3(2.f, 2.f, 2.f);
		pca.addVector(vec);
		direction = pca.getPrincipalDirection();
		EXPECT_EQ(glm::abs(direction), norm);

		vec = glm::vec3(3.f, 3.f, 3.f);
		pca.addVector(vec);
		direction = pca.getPrincipalDirection();
		EXPECT_EQ(glm::abs(direction), norm);

		vec = glm::vec3(4.f, 4.f, 4.f);
		pca.addVector(vec);
		direction = pca.getPrincipalDirection();
		EXPECT_EQ(glm::abs(direction), norm);
	}

	TEST(Direction, OppositeDirections) {
		OnlinePCA pca;

		glm::vec3 dir1(1.f, 1.f, 1.f);
		pca.addVector(dir1);
		glm::vec3 direction = pca.getPrincipalDirection();
		EXPECT_EQ(direction, glm::normalize(dir1));

		glm::vec3 dir2(-1.f, -1.f, -1.f);
		pca.addVector(dir2);
		direction = pca.getPrincipalDirection();
		EXPECT_TRUE(direction == glm::normalize(dir1) || direction == glm::normalize(dir2));
	}

	TEST(Direction, MultipleClearDirections) {
		OnlinePCA pca;

		glm::vec3 dir1(1.f, 0.f, 0.f);
		pca.addVector(dir1);
		glm::vec3 direction = pca.getPrincipalDirection();
		EXPECT_EQ(direction, dir1);

		glm::vec3 dir2 = glm::vec3(0.f, 1.f, 0.f);

		for (int i = 0; i < 100; i++) {
			pca.addVector(dir2);
			pca.addVector(dir1);

			direction = pca.getPrincipalDirection();
			EXPECT_TRUE(direction == dir1 || direction == dir2);
		}
	}
};
