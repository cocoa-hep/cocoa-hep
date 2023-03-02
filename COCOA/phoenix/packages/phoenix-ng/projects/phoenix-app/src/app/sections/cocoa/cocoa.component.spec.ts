import { ComponentFixture, TestBed } from '@angular/core/testing';

import { ScdComponent } from './cocoa.component';

describe('ScdComponent', () => {
  let component: ScdComponent;
  let fixture: ComponentFixture<ScdComponent>;

  beforeEach(async () => {
    await TestBed.configureTestingModule({
      declarations: [ScdComponent],
    }).compileComponents();
  });

  beforeEach(() => {
    fixture = TestBed.createComponent(ScdComponent);
    component = fixture.componentInstance;
    fixture.detectChanges();
  });

  it('should create', () => {
    expect(component).toBeTruthy();
  });
});
